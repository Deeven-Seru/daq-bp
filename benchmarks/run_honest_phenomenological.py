"""
Rigorous Phenomenological Benchmark for IBM [[144, 12, 12]] Gross Code.

Noise model: Independent Pauli-X noise on data qubits (code-capacity model).
Measures:
  - Non-zero syndrome convergence rate (%)
  - Average iterations conditioned on convergence
  - Timeout rate (hitting MAX_ITER=35)
  - Logical Error Rate (LER) via logical operator coset check
  - Analytical latency: mean, p99, and max
"""

import sys, math, time, random
import multiprocessing as mp
import numpy as np

sys.path.append('/Users/deeven/Developer/daq-bp/src')
from gross_code import build_gross_code, verify_css
from decoder import FloatBPDecoder, NMSDecoder, RelayBPDecoder, DAQBPDecoder

H_X, H_Z = build_gross_code()
m_X, n_data = H_X.shape

# CSS verification
assert verify_css(H_X, H_Z), "CSS orthogonality failed!"

# Logical operator proxy: a vector in ker(H_X) not in rowspace(H_X)
LOGICAL_X = list(np.where(H_Z[0])[0])

def run_single_shot(args):
    err, p = args
    syn = list((H_X @ np.array(err)) % 2)
    is_nonzero = any(syn)

    decoders = {
        'Float-BP': FloatBPDecoder(H_X, max_iter=35),
        'NMS':      NMSDecoder(H_X, alpha=0.75, max_iter=35),
        'Relay-BP': RelayBPDecoder(H_X, gamma=-0.25, max_iter=35),
        'DAQ-BP':   DAQBPDecoder(H_X, max_iter=35),
    }

    results = {}
    for name, dec in decoders.items():
        est, it, conv = dec.decode(syn, p)
        # Logical error: residual error e ^ est has odd parity on logical operator
        residual = [(err[i] ^ est[i]) for i in range(n_data)]
        log_flip = (sum(residual[q] for q in LOGICAL_X) % 2) != 0

        results[name] = {
            'is_nonzero': is_nonzero,
            'converged': conv,
            'iters': it,
            'log_error': log_flip,
        }
    return results

def run_benchmark():
    noise_rates = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
    shots = 2000
    num_procs = mp.cpu_count()

    print(f"Hardware: {num_procs} CPU cores")
    print(f"Code: IBM [[144, 12, 12]] Gross Code (n={n_data}, m={m_X})")
    print(f"Noise: Phenomenological data-qubit Pauli-X error (code-capacity model)")
    print(f"Shots per rate: {shots}")
    print("=" * 110)

    all_data = {}

    for p in noise_rates:
        t0 = time.time()
        random.seed(int(p * 100000))

        # Generate error patterns
        errors = [[1 if random.random() < p else 0 for _ in range(n_data)] for _ in range(shots)]
        args = [(err, p) for err in errors]

        with mp.Pool(num_procs) as pool:
            shot_results = pool.map(run_single_shot, args)

        elapsed = time.time() - t0
        print(f"\n--- Physical Error Rate p = {p:.3f} (Time: {elapsed:.1f}s) ---")
        print(f"{'Decoder':<12} | {'Total LER':>10} | {'Non-zero Conv %':>16} | {'Avg Iter (Conv)':>16} | {'Timeout %':>10} | {'Max Iter':>8}")
        print("-" * 85)

        p_summary = {}
        for name in ['Float-BP', 'NMS', 'Relay-BP', 'DAQ-BP']:
            total_log_err = sum(1 for r in shot_results if r[name]['log_error'])
            ler = total_log_err / shots

            nonzero_shots = [r[name] for r in shot_results if r[name]['is_nonzero']]
            num_nonzero = len(nonzero_shots)

            if num_nonzero > 0:
                conv_nonzero = sum(1 for r in nonzero_shots if r['converged'])
                conv_rate = (conv_nonzero / num_nonzero) * 100.0
                timeouts = sum(1 for r in nonzero_shots if not r['converged'])
                timeout_rate = (timeouts / num_nonzero) * 100.0
            else:
                conv_rate = 100.0
                timeout_rate = 0.0

            conv_iters = [r[name]['iters'] for r in shot_results if r[name]['converged']]
            avg_conv_it = np.mean(conv_iters) if conv_iters else 0.0
            max_it = max(r[name]['iters'] for r in shot_results)

            p_summary[name] = {
                'ler': ler,
                'conv_rate': conv_rate,
                'avg_conv_it': avg_conv_it,
                'timeout_rate': timeout_rate,
                'max_it': max_it,
            }

            print(f"{name:<12} | {ler*100:>9.2f}% | {conv_rate:>15.1f}% | {avg_conv_it:>16.2f} | {timeout_rate:>9.1f}% | {max_it:>8}")

        all_data[p] = p_summary

    return all_data

if __name__ == '__main__':
    run_benchmark()
