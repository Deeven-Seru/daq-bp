"""
Decoder implementations for quantum LDPC codes:
- DAQ-BP (Degeneracy-Aware Quantized 5-bit Integer BP)
- Float-BP (Standard Sum-Product Floating-Point BP)
- NMS (Normalized Min-Sum)
- Relay-BP (Float-BP with memory damping)
"""

import math
from typing import List, Tuple
import numpy as np


class TannerGraph:
    """Represents Tanner graph adjacency structures."""

    def __init__(self, H: np.ndarray):
        self.H = H
        self.num_checks, self.num_bits = H.shape
        self.check_adj = [list(np.where(H[c])[0]) for c in range(self.num_checks)]
        self.bit_adj = [list(np.where(H[:, v])[0]) for v in range(self.num_bits)]
        self.vc_idx = {}
        for v in range(self.num_bits):
            for idx, c in enumerate(self.bit_adj[v]):
                self.vc_idx[(v, c)] = idx


class DAQBPDecoder:
    """
    5-bit Integer Quantized Belief Propagation with:
      - Iteration-Adaptive Degeneracy Damping: alpha_deg(k) = >> (k // 10)
      - Anti-Trapping Momentum Operator: sign-inversion detection & halving
    """

    SCALE_Q = 4.0

    def __init__(self, H: np.ndarray, max_iter: int = 35):
        self.graph = TannerGraph(H)
        self.max_iter = max_iter

        # Precompute 5-bit signed integer LUTs
        self.tanh_lut = np.array(
            [round(math.tanh(q / (self.SCALE_Q * 2)) * 32) for q in range(-15, 16)],
            dtype=np.int32,
        )
        self.atanh_lut = np.array(
            [
                round(math.atanh(max(-0.9999, min(0.9999, t / 32.0))) * self.SCALE_Q * 2)
                for t in range(-32, 33)
            ],
            dtype=np.int32,
        )

    def decode(self, syndrome: List[int], p: float) -> Tuple[List[int], int, bool]:
        """Decode binary syndrome. Returns (correction, iterations, converged)."""
        m_X = self.graph.num_checks
        n_data = self.graph.num_bits
        ch_llr = math.log((1.0 - p) / max(p, 1e-12))
        ch_llr_q = int(np.clip(round(ch_llr * self.SCALE_Q), -15, 15))

        v2c = [[ch_llr_q] * len(self.graph.bit_adj[v]) for v in range(n_data)]
        c2v = [[0] * len(self.graph.check_adj[c]) for c in range(m_X)]

        prev_est = [0] * n_data
        trap_cnt = [0] * n_data

        for it in range(1, self.max_iter + 1):
            shift = (it - 1) // 10

            # Check node updates (Integer LUT-based tanh product)
            for c in range(m_X):
                s = syndrome[c]
                vars_c = self.graph.check_adj[c]
                tanh_v = [
                    int(self.tanh_lut[max(-15, min(15, v2c[v][self.graph.vc_idx[(v, c)]])) + 15])
                    for v in vars_c
                ]
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
                    raw_q = int(self.atanh_lut[max(-32, min(32, local_sign * prod)) + 32])
                    damped = raw_q >> shift
                    c2v[c][idx_v] = max(-24, min(24, damped))

            # Variable node updates + Anti-trapping momentum
            est = [0] * n_data
            for v in range(n_data):
                tot = ch_llr_q + sum(c2v[c][self.graph.vc_idx[(v, c)]] for c in self.graph.bit_adj[v])
                bit = 1 if tot < 0 else 0
                est[v] = bit
                if bit != prev_est[v]:
                    trap_cnt[v] += 1
                else:
                    trap_cnt[v] = max(0, trap_cnt[v] - 1)

                for idx_c, c in enumerate(self.graph.bit_adj[v]):
                    raw = ch_llr_q + sum(
                        c2v[cc][self.graph.vc_idx[(v, cc)]]
                        for cc in self.graph.bit_adj[v]
                        if cc != c
                    )
                    if trap_cnt[v] >= 3:
                        raw = raw >> 1
                    v2c[v][idx_c] = max(-24, min(24, raw))

            prev_est = est[:]

            # Check syndrome satisfaction
            syn_check = [0] * m_X
            for c in range(m_X):
                for v in self.graph.check_adj[c]:
                    syn_check[c] ^= est[v]
            if syn_check == syndrome:
                return est, it, True

        return est, self.max_iter, False


class FloatBPDecoder:
    """Standard Sum-Product Belief Propagation (64-bit float)."""

    def __init__(self, H: np.ndarray, max_iter: int = 35):
        self.graph = TannerGraph(H)
        self.max_iter = max_iter

    def decode(self, syndrome: List[int], p: float) -> Tuple[List[int], int, bool]:
        m_X = self.graph.num_checks
        n_data = self.graph.num_bits
        ch_llr = math.log((1.0 - p) / max(p, 1e-12))
        clip_f = 20.0

        v2c = [[ch_llr] * len(self.graph.bit_adj[v]) for v in range(n_data)]
        c2v = [[0.0] * len(self.graph.check_adj[c]) for c in range(m_X)]

        for it in range(1, self.max_iter + 1):
            for c in range(m_X):
                s = syndrome[c]
                vars_c = self.graph.check_adj[c]
                msgs = [v2c[v][self.graph.vc_idx[(v, c)]] for v in vars_c]
                half = [max(-clip_f, min(clip_f, x / 2.0)) for x in msgs]
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
                tot = ch_llr + sum(c2v[c][self.graph.vc_idx[(v, c)]] for c in self.graph.bit_adj[v])
                est[v] = 1 if tot < 0.0 else 0
                for idx_c, c in enumerate(self.graph.bit_adj[v]):
                    v2c[v][idx_c] = ch_llr + sum(
                        c2v[cc][self.graph.vc_idx[(v, cc)]]
                        for cc in self.graph.bit_adj[v]
                        if cc != c
                    )

            syn_check = [0] * m_X
            for c in range(m_X):
                for v in self.graph.check_adj[c]:
                    syn_check[c] ^= est[v]
            if syn_check == syndrome:
                return est, it, True

        return est, self.max_iter, False


class NMSDecoder:
    """Normalized Min-Sum Decoder."""

    def __init__(self, H: np.ndarray, alpha: float = 0.75, max_iter: int = 35):
        self.graph = TannerGraph(H)
        self.alpha = alpha
        self.max_iter = max_iter

    def decode(self, syndrome: List[int], p: float) -> Tuple[List[int], int, bool]:
        m_X = self.graph.num_checks
        n_data = self.graph.num_bits
        ch_llr = math.log((1.0 - p) / max(p, 1e-12))

        v2c = [[ch_llr] * len(self.graph.bit_adj[v]) for v in range(n_data)]
        c2v = [[0.0] * len(self.graph.check_adj[c]) for c in range(m_X)]

        for it in range(1, self.max_iter + 1):
            for c in range(m_X):
                s = syndrome[c]
                vars_c = self.graph.check_adj[c]
                msgs = [v2c[v][self.graph.vc_idx[(v, c)]] for v in vars_c]
                sign_prod = (-1) ** s
                for x in msgs:
                    sign_prod *= (1 if x >= 0 else -1)
                abs_msgs = [abs(x) for x in msgs]
                for idx_v in range(len(vars_c)):
                    others = [abs_msgs[i] for i in range(len(vars_c)) if i != idx_v]
                    min_v = min(others) if others else 0.0
                    local_sign = sign_prod * (1 if msgs[idx_v] >= 0 else -1)
                    c2v[c][idx_v] = self.alpha * local_sign * min_v

            est = [0] * n_data
            for v in range(n_data):
                tot = ch_llr + sum(c2v[c][self.graph.vc_idx[(v, c)]] for c in self.graph.bit_adj[v])
                est[v] = 1 if tot < 0.0 else 0
                for idx_c, c in enumerate(self.graph.bit_adj[v]):
                    v2c[v][idx_c] = ch_llr + sum(
                        c2v[cc][self.graph.vc_idx[(v, cc)]]
                        for cc in self.graph.bit_adj[v]
                        if cc != c
                    )

            syn_check = [0] * m_X
            for c in range(m_X):
                for v in self.graph.check_adj[c]:
                    syn_check[c] ^= est[v]
            if syn_check == syndrome:
                return est, it, True

        return est, self.max_iter, False


class RelayBPDecoder:
    """Relay Belief Propagation Decoder with memory damping."""

    def __init__(self, H: np.ndarray, gamma: float = -0.25, max_iter: int = 35):
        self.graph = TannerGraph(H)
        self.gamma = gamma
        self.max_iter = max_iter

    def decode(self, syndrome: List[int], p: float) -> Tuple[List[int], int, bool]:
        m_X = self.graph.num_checks
        n_data = self.graph.num_bits
        ch_llr = math.log((1.0 - p) / max(p, 1e-12))
        clip_f = 20.0

        v2c = [[ch_llr] * len(self.graph.bit_adj[v]) for v in range(n_data)]
        prev_v2c = [[ch_llr] * len(self.graph.bit_adj[v]) for v in range(n_data)]
        c2v = [[0.0] * len(self.graph.check_adj[c]) for c in range(m_X)]

        for it in range(1, self.max_iter + 1):
            for c in range(m_X):
                s = syndrome[c]
                vars_c = self.graph.check_adj[c]
                msgs = [v2c[v][self.graph.vc_idx[(v, c)]] for v in vars_c]
                half = [max(-clip_f, min(clip_f, x / 2.0)) for x in msgs]
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
                tot = ch_llr + sum(c2v[c][self.graph.vc_idx[(v, c)]] for c in self.graph.bit_adj[v])
                est[v] = 1 if tot < 0.0 else 0
                for idx_c, c in enumerate(self.graph.bit_adj[v]):
                    new_val = ch_llr + sum(
                        c2v[cc][self.graph.vc_idx[(v, cc)]]
                        for cc in self.graph.bit_adj[v]
                        if cc != c
                    )
                    blended = (1.0 - self.gamma) * new_val + self.gamma * prev_v2c[v][idx_c]
                    prev_v2c[v][idx_c] = v2c[v][idx_c]
                    v2c[v][idx_c] = blended

            syn_check = [0] * m_X
            for c in range(m_X):
                for v in self.graph.check_adj[c]:
                    syn_check[c] ^= est[v]
            if syn_check == syndrome:
                return est, it, True

        return est, self.max_iter, False
