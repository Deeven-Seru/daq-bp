"""
Comprehensive Benchmark on IBM [[144, 12, 12]] Gross Code.
Compares:
  1. Float-BP (Standalone)
  2. NMS (Standalone)
  3. Relay-BP (Standalone)
  4. DAQ-BP (Standalone 5-bit Integer)
  5. Standard BP-OSD (Float MS + OSD-CS)
  6. DAQ-BP + OSD Fallback (5-bit Integer Front-End + OSD Fallback)
"""

import sys, math, time, random
import multiprocessing as mp
import numpy as np

sys.path.append('/Users/deeven/Developer/daq-bp/src')
from gross_code import build_gross_code, verify_css
from decoder import FloatBPDecoder, NMSDecoder, RelayBPDecoder, DAQBPDecoder
from ldpc import bposd_decoder

H_X, H_Z = build_gross_code()
m_X, n_data = H_X.shape
assert verify_css(H_X, H_Z), "CSS check failed!"
LOGICAL_X = list(np.where(H_Z[0])[0])

def run_single_shot(args):
    err, p = args
    syn = (H_X @ np.array(err)) % 2
    syn_list = list(syn)
    is_nonzero = bool(np.any(syn))

    # Standalone decoders
    float_dec = FloatBPDecoder(H_X, max_iter=35)
    nms_dec   = NMSDecoder(H_X, alpha=0.75, max_iter=35)
    relay_dec = RelayBPDecoder(H_X, gamma=-0.25, max_iter=35)
    daq_dec   = DAQBPDecoder(H_X, max_iter=35)

    f_est, f_it, f_conv = float_dec.decode(syn_list, p)
    n_est, n_it, n_conv = nms_dec.decode(syn_list, p)
    r_est, r_it, r_conv = relay_dec.decode(syn_list, p)
    d_est, d_it, d_conv = daq_dec.decode(syn_list, p)

    # Standard BP-OSD
    bposd = bposd_decoder(
        H_X,
        error_rate=max(p, 1e-4),
        max_iter=35,
        bp_method='ms',
        ms_scaling_factor=0.75,
        osd_method='osd_cs',
        osd_order=10
    )
    bposd_est = bposd.decode(syn)

    # Hybrid DAQ-BP + OSD Fallback:
    # If DAQ-BP converges, use its estimate; otherwise fall back to OSD
    if d_conv:
        hybrid_est = d_est
        used_osd = False
    else:
        hybrid_est = bposd_est
        used_osd = True

    def check_log_error(est):
        res = [(err[i] ^ int(est[i])) for i in range(n_data)]
        return (sum(res[q] for q in LOGICAL_X) % 2) != 0

    return {
        'is_nonzero': is_nonzero,
        'Float-BP':    {'log_err': check_log_error(f_est), 'conv': f_conv, 'iter': f_it},
        'NMS':         {'log_err': check_log_error(n_est), 'conv': n_conv, 'iter': n_it},
        'Relay-BP':    {'log_err': check_log_error(r_est), 'conv': r_conv, 'iter': r_it},
        'DAQ-BP':      {'log_err': check_log_error(d_est), 'conv': d_conv, 'iter': d_it},
        'BP-OSD':      {'log_err': check_log_error(bposd_est), 'conv': True, 'iter': 35},
        'DAQ-BP+OSD':  {'log_err': check_log_error(hybrid_est), 'conv': True, 'used_osd': used_osd, 'iter': d_it if not used_osd else 35},
    }

def run_benchmark():
    noise_rates = [0.005, 0.010, 0.020, 0.030, 0.040, 0.050]
    shots = 1000
    num_procs = mp.cpu_count()

    print(f"Executing Multi-Decoder Benchmark with Standard BP-OSD Reference ({num_procs} cores)")
    print(f"IBM [[144, 12, 12]] Gross Code (n={n_data}, m={m_X})")
    print(f"Monte Carlo Shots per Rate: {shots}")
    print("=" * 115)

    summary = {}

    for p in noise_rates:
        t0 = time.time()
        random.seed(int(p * 100000))
        errors = [[1 if random.random() < p else 0 for _ in range(n_data)] for _ in range(shots)]
        args = [(err, p) for err in errors]

        with mp.Pool(num_procs) as pool:
            results = pool.map(run_single_shot, args)

        elapsed = time.time() - t0
        print(f"\n--- Physical Error Rate p = {p:.3f} (Time: {elapsed:.1f}s) ---")
        print(f"{'Decoder':<15} | {'Logical Error Rate':>20} | {'Non-Zero Conv %':>18} | {'Avg Converged Iter':>20} | {'OSD Fallback %':>16}")
        print("-" * 98)

        p_res = {}
        for name in ['Float-BP', 'NMS', 'Relay-BP', 'DAQ-BP', 'BP-OSD', 'DAQ-BP+OSD']:
            total_log_err = sum(1 for r in results if r[name]['log_err'])
            ler = total_log_err / shots

            nonzero = [r[name] for r in results if r['is_nonzero']]
            if nonzero:
                conv_rate = (sum(1 for r in nonzero if r['conv']) / len(nonzero)) * 100.0
            else:
                conv_rate = 100.0

            conv_iters = [r[name]['iter'] for r in results if r[name]['conv']]
            avg_it = np.mean(conv_iters) if conv_iters else 0.0

            if name == 'DAQ-BP+OSD':
                osd_fallback_pct = (sum(1 for r in results if r[name]['used_osd']) / shots) * 100.0
            else:
                osd_fallback_pct = 100.0 if name == 'BP-OSD' else 0.0

            p_res[name] = {
                'ler': ler,
                'conv_rate': conv_rate,
                'avg_it': avg_it,
                'osd_fallback_pct': osd_fallback_pct
            }

            osd_str = f"{osd_fallback_pct:>15.1f}%" if name in ['BP-OSD', 'DAQ-BP+OSD'] else "N/A"
            print(f"{name:<15} | {ler*100:>19.2f}% | {conv_rate:>17.1f}% | {avg_it:>20.2f} | {osd_str:>16}")

        summary[p] = p_res

    return summary

if __name__ == '__main__':
    run_benchmark()
