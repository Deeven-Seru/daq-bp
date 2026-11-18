"""
Parallel circuit-level benchmark for IBM [[144, 12, 12]] Gross code using Stim.
Evaluates 4 decoders:
  1. Float-BP
  2. NMS (Normalized Min-Sum)
  3. Relay-BP
  4. DAQ-BP (Integer-only 5-bit, our proposal)

Uses multiprocessing.Pool for 8x speedup.
Measures true logical error rate (observable parity) + average iterations.
"""

import stim
import math
import time
import sys
import multiprocessing as mp
import numpy as np

# ---------------------------------------------------------------------------
# 1. Gross code construction & CSS check
# ---------------------------------------------------------------------------

def build_gross_code():
    l, m = 12, 6
    L = l * m  # 72

    def index(xi, yj):
        return (xi % l) * m + (yj % m)

    A = np.zeros((L, L), dtype=np.int8)
    B = np.zeros((L, L), dtype=np.int8)
    for i in range(l):
        for j in range(m):
            row = index(i, j)
            A[row, index(i+3, j  )] = 1
            A[row, index(i,   j+1)] = 1
            A[row, index(i,   j+2)] = 1
            B[row, index(i,   j+3)] = 1
            B[row, index(i+1, j  )] = 1
            B[row, index(i+2, j  )] = 1

    H_X = np.concatenate([A, B], axis=1)
    H_Z = np.concatenate([B.T, A.T], axis=1)
    return H_X, H_Z

# Precompute adjacency
H_X, H_Z = build_gross_code()
m_X, n_data = H_X.shape
check_adj = [list(np.where(H_X[c])[0]) for c in range(m_X)]
var_adj   = [list(np.where(H_X[:, v])[0]) for v in range(n_data)]
vc_idx = {}
for v in range(n_data):
    for idx, c in enumerate(var_adj[v]):
        vc_idx[(v, c)] = idx

# Logical X proxy support (from ker(H_X))
LOGICAL_X_SUPPORT = list(np.where(H_Z[0])[0])

# Precompute LUTs for DAQ-BP
SCALE_Q = 4.0
TANH_LUT = np.array([round(math.tanh(q / (SCALE_Q * 2)) * 32) for q in range(-15, 16)], dtype=np.int32)
ATANH_LUT = np.array([round(math.atanh(max(-0.9999, min(0.9999, t / 32.0))) * SCALE_Q * 2) for t in range(-32, 33)], dtype=np.int32)

# ---------------------------------------------------------------------------
# 2. Stim circuit builder
# ---------------------------------------------------------------------------

def build_circuit(H_X, rounds, p, obs_support):
    m_X, n = H_X.shape
    anc_off = n
    c = stim.Circuit()
    c.append("R", list(range(anc_off, anc_off + m_X)))

    for rnd in range(rounds):
        if rnd > 0:
            c.append("R", list(range(anc_off, anc_off + m_X)))
        for ci in range(m_X):
            anc = anc_off + ci
            for dq in check_adj[ci]:
                c.append("CNOT", [anc, int(dq)])
                c.append("DEPOLARIZE2", [anc, int(dq)], max(1e-9, p / 10.0))
        c.append("DEPOLARIZE1", list(range(n)), max(1e-9, p))
        for ci in range(m_X):
            c.append("X_ERROR", [anc_off + ci], max(1e-9, p))
        c.append("M", list(range(anc_off, anc_off + m_X)))
        if rnd == 0:
            for ci in range(m_X):
                c.append("DETECTOR", [stim.target_rec(-(m_X - ci))])
        else:
            for ci in range(m_X):
                c.append("DETECTOR", [
                    stim.target_rec(-(m_X - ci)),
                    stim.target_rec(-(m_X - ci) - m_X),
                ])

    c.append("M", list(range(n)))
    obs = [stim.target_rec(-(n - int(q))) for q in obs_support]
    c.append("OBSERVABLE_INCLUDE", obs, 0)
    return c

# ---------------------------------------------------------------------------
# 3. Decoder implementations
# ---------------------------------------------------------------------------

MAX_ITER = 35
CLIP_F = 20.0

def decode_float_bp(syndrome, ch_llr):
    v2c = [[ch_llr] * len(var_adj[v]) for v in range(n_data)]
    c2v = [[0.0] * len(check_adj[c]) for c in range(m_X)]

    for it in range(1, MAX_ITER + 1):
        for c in range(m_X):
            s = syndrome[c]
            vars_c = check_adj[c]
            msgs = [v2c[v][vc_idx[(v, c)]] for v in vars_c]
            half = [max(-CLIP_F, min(CLIP_F, x / 2.0)) for x in msgs]
            tanh_v = [max(-1 + 1e-9, min(1 - 1e-9, math.tanh(x))) for x in half]
            prod_all = 1.0
            for t in tanh_v:
                prod_all *= t
            prod_all *= (-1) ** s
            for idx_v in range(len(vars_c)):
                ex = prod_all / tanh_v[idx_v]
                ex = max(-1 + 1e-9, min(1 - 1e-9, ex))
                c2v[c][idx_v] = 2.0 * math.atanh(ex)

        est = [0] * n_data
        for v in range(n_data):
            tot = ch_llr + sum(c2v[c][vc_idx[(v, c)]] for c in var_adj[v])
            est[v] = 1 if tot < 0.0 else 0
            for idx_c, c in enumerate(var_adj[v]):
                v2c[v][idx_c] = ch_llr + sum(c2v[cc][vc_idx[(v, cc)]] for cc in var_adj[v] if cc != c)

        syn_check = [0] * m_X
        for c in range(m_X):
            for v in check_adj[c]:
                syn_check[c] ^= est[v]
        if syn_check == syndrome:
            return est, it, True

    return est, MAX_ITER, False


def decode_nms(syndrome, ch_llr, alpha=0.75):
    v2c = [[ch_llr] * len(var_adj[v]) for v in range(n_data)]
    c2v = [[0.0] * len(check_adj[c]) for c in range(m_X)]

    for it in range(1, MAX_ITER + 1):
        for c in range(m_X):
            s = syndrome[c]
            vars_c = check_adj[c]
            msgs = [v2c[v][vc_idx[(v, c)]] for v in vars_c]
            sign_prod = (-1) ** s
            for x in msgs:
                sign_prod *= (1 if x >= 0 else -1)
            abs_msgs = [abs(x) for x in msgs]
            for idx_v in range(len(vars_c)):
                others = [abs_msgs[i] for i in range(len(vars_c)) if i != idx_v]
                min_v = min(others) if others else 0.0
                local_sign = sign_prod * (1 if msgs[idx_v] >= 0 else -1)
                c2v[c][idx_v] = alpha * local_sign * min_v

        est = [0] * n_data
        for v in range(n_data):
            tot = ch_llr + sum(c2v[c][vc_idx[(v, c)]] for c in var_adj[v])
            est[v] = 1 if tot < 0.0 else 0
            for idx_c, c in enumerate(var_adj[v]):
                v2c[v][idx_c] = ch_llr + sum(c2v[cc][vc_idx[(v, cc)]] for cc in var_adj[v] if cc != c)

        syn_check = [0] * m_X
        for c in range(m_X):
            for v in check_adj[c]:
                syn_check[c] ^= est[v]
        if syn_check == syndrome:
            return est, it, True

    return est, MAX_ITER, False


def decode_relay_bp(syndrome, ch_llr, gamma=-0.25):
    v2c = [[ch_llr] * len(var_adj[v]) for v in range(n_data)]
    prev_v2c = [[ch_llr] * len(var_adj[v]) for v in range(n_data)]
    c2v = [[0.0] * len(check_adj[c]) for c in range(m_X)]

    for it in range(1, MAX_ITER + 1):
        for c in range(m_X):
            s = syndrome[c]
            vars_c = check_adj[c]
            msgs = [v2c[v][vc_idx[(v, c)]] for v in vars_c]
            half = [max(-CLIP_F, min(CLIP_F, x / 2.0)) for x in msgs]
            tanh_v = [max(-1 + 1e-9, min(1 - 1e-9, math.tanh(x))) for x in half]
            prod_all = 1.0
            for t in tanh_v:
                prod_all *= t
            prod_all *= (-1) ** s
            for idx_v in range(len(vars_c)):
                ex = prod_all / tanh_v[idx_v]
                ex = max(-1 + 1e-9, min(1 - 1e-9, ex))
                c2v[c][idx_v] = 2.0 * math.atanh(ex)

        est = [0] * n_data
        for v in range(n_data):
            tot = ch_llr + sum(c2v[c][vc_idx[(v, c)]] for c in var_adj[v])
            est[v] = 1 if tot < 0.0 else 0
            for idx_c, c in enumerate(var_adj[v]):
                new_val = ch_llr + sum(c2v[cc][vc_idx[(v, cc)]] for cc in var_adj[v] if cc != c)
                blended = (1.0 - gamma) * new_val + gamma * prev_v2c[v][idx_c]
                prev_v2c[v][idx_c] = v2c[v][idx_c]
                v2c[v][idx_c] = blended

        syn_check = [0] * m_X
        for c in range(m_X):
            for v in check_adj[c]:
                syn_check[c] ^= est[v]
        if syn_check == syndrome:
            return est, it, True

    return est, MAX_ITER, False


def decode_daq_bp(syndrome, ch_llr):
    ch_llr_q = int(np.clip(round(ch_llr * SCALE_Q), -15, 15))
    v2c = [[ch_llr_q] * len(var_adj[v]) for v in range(n_data)]
    c2v = [[0] * len(check_adj[c]) for c in range(m_X)]

    prev_est = [0] * n_data
    trap_cnt = [0] * n_data

    for it in range(1, MAX_ITER + 1):
        shift = (it - 1) // 10

        for c in range(m_X):
            s = syndrome[c]
            vars_c = check_adj[c]
            tanh_v = [int(TANH_LUT[max(-15, min(15, v2c[v][vc_idx[(v, c)]])) + 15]) for v in vars_c]
            signs = [1 if t >= 0 else -1 for t in tanh_v]
            abs_t = [abs(t) for t in tanh_v]
            prod_sign = (-1) ** s
            for s_ in signs:
                prod_sign *= s_
            for idx_v in range(len(vars_c)):
                prod = 32
                for idx2 in range(len(vars_c)):
                    if idx2 != idx_v:
                        prod = (prod * abs_t[idx2]) >> 5
                        prod = max(0, min(31, prod))
                local_sign = prod_sign * signs[idx_v]
                raw_q = int(ATANH_LUT[max(-32, min(32, local_sign * prod)) + 32])
                damped = raw_q >> shift
                c2v[c][idx_v] = max(-24, min(24, damped))

        est = [0] * n_data
        for v in range(n_data):
            tot = ch_llr_q + sum(c2v[c][vc_idx[(v, c)]] for c in var_adj[v])
            bit = 1 if tot < 0 else 0
            est[v] = bit
            if bit != prev_est[v]:
                trap_cnt[v] += 1
            else:
                trap_cnt[v] = max(0, trap_cnt[v] - 1)
            for idx_c, c in enumerate(var_adj[v]):
                raw = ch_llr_q + sum(c2v[cc][vc_idx[(v, cc)]] for cc in var_adj[v] if cc != c)
                if trap_cnt[v] >= 3:
                    raw = raw >> 1
                v2c[v][idx_c] = max(-24, min(24, raw))

        prev_est = est[:]

        syn_check = [0] * m_X
        for c in range(m_X):
            for v in check_adj[c]:
                syn_check[c] ^= est[v]
        if syn_check == syndrome:
            return est, it, True

    return est, MAX_ITER, False

# ---------------------------------------------------------------------------
# 4. Worker function for parallel shot processing
# ---------------------------------------------------------------------------

def worker_process_chunk(args):
    """Processes a chunk of (syndrome, actual_logical_flip, ch_llr)."""
    shots_chunk, ch_llr = args
    results = {
        'Float-BP': {'correct': 0, 'iters': 0, 'converged': 0},
        'NMS':      {'correct': 0, 'iters': 0, 'converged': 0},
        'Relay-BP': {'correct': 0, 'iters': 0, 'converged': 0},
        'DAQ-BP':   {'correct': 0, 'iters': 0, 'converged': 0},
    }

    decoders = {
        'Float-BP': decode_float_bp,
        'NMS':      decode_nms,
        'Relay-BP': decode_relay_bp,
        'DAQ-BP':   decode_daq_bp,
    }

    for syn, actual_flip in shots_chunk:
        for name, dec_fn in decoders.items():
            est, iters, conv = dec_fn(syn, ch_llr)
            results[name]['iters'] += iters
            if conv:
                results[name]['converged'] += 1
            corr_parity = sum(est[q] for q in LOGICAL_X_SUPPORT) % 2
            if corr_parity == int(actual_flip):
                results[name]['correct'] += 1

    return results

# ---------------------------------------------------------------------------
# 5. Main benchmark harness
# ---------------------------------------------------------------------------

def run_multiprocess_benchmark():
    num_procs = mp.cpu_count()
    print(f"Hardware: {num_procs} CPU cores detected.")
    print("Code: IBM [[144, 12, 12]] Gross code")
    print(f"CSS Orthogonality: {not np.any((H_X @ H_Z.T) % 2)}")
    print(f"Decoders: Float-BP, NMS (alpha=0.75), Relay-BP (gamma=-0.25), DAQ-BP (5-bit integer)")

    noise_rates = [0.001, 0.002, 0.004, 0.006, 0.008, 0.010]
    shots = 1000
    rounds = 3

    print(f"Evaluation: {rounds} syndrome extraction rounds (Stim circuit-level noise), {shots} shots per rate.")
    print("-" * 88)
    hdr = f"{'p':>6} | {'Float-BP LER':>14} | {'NMS LER':>14} | {'Relay-BP LER':>14} | {'DAQ-BP (Ours)':>14} | {'Time':>6}"
    print(hdr)
    print("-" * 88)

    summary_table = []

    pool = mp.Pool(processes=num_procs)

    for p in noise_rates:
        t0 = time.time()
        circ = build_circuit(H_X, rounds, p, LOGICAL_X_SUPPORT)
        det_sampler = circ.compile_detector_sampler(seed=1337)
        raw_det, raw_obs = det_sampler.sample(shots=shots, bit_packed=False, separate_observables=True)

        # Build syndromes
        syndromes = np.zeros((shots, m_X), dtype=np.int8)
        for r in range(rounds):
            syndromes ^= raw_det[:, r * m_X:(r+1) * m_X].astype(np.int8)

        shot_pairs = [(list(syndromes[i]), bool(raw_obs[i, 0])) for i in range(shots)]
        ch_llr = math.log((1.0 - p) / max(p, 1e-12))

        # Chunk for multiprocessing
        chunk_size = math.ceil(shots / num_procs)
        chunks = [(shot_pairs[i:i + chunk_size], ch_llr) for i in range(0, shots, chunk_size)]

        chunk_results = pool.map(worker_process_chunk, chunks)

        # Aggregate results
        agg = {
            'Float-BP': {'correct': 0, 'iters': 0, 'converged': 0},
            'NMS':      {'correct': 0, 'iters': 0, 'converged': 0},
            'Relay-BP': {'correct': 0, 'iters': 0, 'converged': 0},
            'DAQ-BP':   {'correct': 0, 'iters': 0, 'converged': 0},
        }
        for cr in chunk_results:
            for name in agg:
                agg[name]['correct']   += cr[name]['correct']
                agg[name]['iters']     += cr[name]['iters']
                agg[name]['converged'] += cr[name]['converged']

        elapsed = time.time() - t0

        row_str = f"{p:>6.3f} | "
        p_row = {'p': p, 'elapsed': elapsed}
        for name in ['Float-BP', 'NMS', 'Relay-BP', 'DAQ-BP']:
            ler = 1.0 - (agg[name]['correct'] / shots)
            avg_it = agg[name]['iters'] / shots
            p_row[name + '_ler'] = ler
            p_row[name + '_iters'] = avg_it
            p_row[name + '_conv'] = agg[name]['converged'] / shots
            row_str += f"{ler*100:>10.2f}% ({avg_it:>3.1f}) | "
        row_str += f"{elapsed:>5.1f}s"
        print(row_str)
        sys.stdout.flush()
        summary_table.append(p_row)

    pool.close()
    pool.join()
    print("-" * 88)
    return summary_table

if __name__ == '__main__':
    run_multiprocess_benchmark()
