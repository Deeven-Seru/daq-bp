"""
Comprehensive Benchmark on the authentic IBM [[144, 12, 12]] Gross Code.
Compares:
1. Standard Floating-Point Sum-Product BP
2. Classical Normalized Min-Sum (NMS, alpha=0.75)
3. Relay-BP Baseline (Float, memory perturbation)
4. DAQ-BP (Our Novel 5-bit Integer Quantized Algorithm, 0 DSPs)
"""

import math
import random
import time
from typing import List, Tuple
from gross_code_144 import build_gross_code_matrices

class TannerGraph:
    def __init__(self, H: List[List[int]]):
        self.H = H
        self.num_checks = len(H)
        self.num_bits = len(H[0])
        self.check_adj: List[List[int]] = [[] for _ in range(self.num_checks)]
        self.bit_adj: List[List[int]] = [[] for _ in range(self.num_bits)]
        for c in range(self.num_checks):
            for v in range(self.num_bits):
                if H[c][v] == 1:
                    self.check_adj[c].append(v)
                    self.bit_adj[v].append(c)

class GrossCodeDecoders:
    @staticmethod
    def standard_float_bp(graph: TannerGraph, syndrome: List[int], p: float, max_iter: int = 35) -> Tuple[bool, int]:
        eps = 1e-12
        p_safe = max(eps, min(1.0 - eps, p))
        init_llr = math.log((1.0 - p_safe) / p_safe)
        
        m_c2v = { (c, v): 0.0 for c in range(graph.num_checks) for v in graph.check_adj[c] }
        m_v2c = { (v, c): init_llr for v in range(graph.num_bits) for c in graph.bit_adj[v] }
        
        for it in range(1, max_iter + 1):
            for c in range(graph.num_checks):
                for v in graph.check_adj[c]:
                    prod = 1.0
                    for v_prime in graph.check_adj[c]:
                        if v_prime == v: continue
                        clamped = max(-25.0, min(25.0, m_v2c[(v_prime, c)] / 2.0))
                        prod *= math.tanh(clamped)
                    prod = max(-0.9999999, min(0.9999999, prod))
                    sign = -1.0 if syndrome[c] == 1 else 1.0
                    m_c2v[(c, v)] = 2.0 * math.atanh(sign * prod)
            
            candidate_error = [0] * graph.num_bits
            for v in range(graph.num_bits):
                total_llr = init_llr + sum(m_c2v[(c, v)] for c in graph.bit_adj[v])
                candidate_error[v] = 1 if total_llr < 0.0 else 0
                for c in graph.bit_adj[v]:
                    m_v2c[(v, c)] = init_llr + sum(m_c2v[(c_prime, v)] for c_prime in graph.bit_adj[v] if c_prime != c)
            
            # Check syndrome
            satisfied = True
            for c in range(graph.num_checks):
                if sum(candidate_error[v] for v in graph.check_adj[c]) % 2 != syndrome[c]:
                    satisfied = False
                    break
            if satisfied:
                return True, it
        return False, max_iter

    @staticmethod
    def classical_nms(graph: TannerGraph, syndrome: List[int], p: float, max_iter: int = 35) -> Tuple[bool, int]:
        eps = 1e-12
        p_safe = max(eps, min(1.0 - eps, p))
        init_llr = math.log((1.0 - p_safe) / p_safe)
        alpha = 0.75
        
        m_c2v = { (c, v): 0.0 for c in range(graph.num_checks) for v in graph.check_adj[c] }
        m_v2c = { (v, c): init_llr for v in range(graph.num_bits) for c in graph.bit_adj[v] }
        
        for it in range(1, max_iter + 1):
            for c in range(graph.num_checks):
                for v in graph.check_adj[c]:
                    sgn = -1 if syndrome[c] == 1 else 1
                    min_abs = float('inf')
                    for v_prime in graph.check_adj[c]:
                        if v_prime == v: continue
                        msg = m_v2c[(v_prime, c)]
                        if msg < 0: sgn *= -1
                        if abs(msg) < min_abs: min_abs = abs(msg)
                    m_c2v[(c, v)] = sgn * alpha * min_abs
            
            candidate_error = [0] * graph.num_bits
            for v in range(graph.num_bits):
                total_llr = init_llr + sum(m_c2v[(c, v)] for c in graph.bit_adj[v])
                candidate_error[v] = 1 if total_llr < 0.0 else 0
                for c in graph.bit_adj[v]:
                    m_v2c[(v, c)] = init_llr + sum(m_c2v[(c_prime, v)] for c_prime in graph.bit_adj[v] if c_prime != c)
            
            satisfied = True
            for c in range(graph.num_checks):
                if sum(candidate_error[v] for v in graph.check_adj[c]) % 2 != syndrome[c]:
                    satisfied = False
                    break
            if satisfied:
                return True, it
        return False, max_iter

    @staticmethod
    def relay_bp_float(graph: TannerGraph, syndrome: List[int], p: float, max_iter: int = 35) -> Tuple[bool, int]:
        """Relay-BP Baseline (Float with memory perturbation as in Müller et al. 2025)."""
        eps = 1e-12
        p_safe = max(eps, min(1.0 - eps, p))
        init_llr = math.log((1.0 - p_safe) / p_safe)
        
        m_c2v = { (c, v): 0.0 for c in range(graph.num_checks) for v in graph.check_adj[c] }
        m_v2c = { (v, c): init_llr for v in range(graph.num_bits) for c in graph.bit_adj[v] }
        prev_m_v2c = { (v, c): init_llr for v in range(graph.num_bits) for c in graph.bit_adj[v] }
        
        # Memory damping factor gamma
        gamma = -0.25 # Negative memory perturbation to break symmetries
        
        for it in range(1, max_iter + 1):
            for c in range(graph.num_checks):
                for v in graph.check_adj[c]:
                    prod = 1.0
                    for v_prime in graph.check_adj[c]:
                        if v_prime == v: continue
                        clamped = max(-25.0, min(25.0, m_v2c[(v_prime, c)] / 2.0))
                        prod *= math.tanh(clamped)
                    prod = max(-0.9999999, min(0.9999999, prod))
                    sign = -1.0 if syndrome[c] == 1 else 1.0
                    m_c2v[(c, v)] = 2.0 * math.atanh(sign * prod)
            
            candidate_error = [0] * graph.num_bits
            for v in range(graph.num_bits):
                total_llr = init_llr + sum(m_c2v[(c, v)] for c in graph.bit_adj[v])
                candidate_error[v] = 1 if total_llr < 0.0 else 0
                for c in graph.bit_adj[v]:
                    new_val = init_llr + sum(m_c2v[(c_prime, v)] for c_prime in graph.bit_adj[v] if c_prime != c)
                    # Relay memory update
                    blended = (1.0 - gamma) * new_val + gamma * prev_m_v2c[(v, c)]
                    prev_m_v2c[(v, c)] = m_v2c[(v, c)]
                    m_v2c[(v, c)] = blended
            
            satisfied = True
            for c in range(graph.num_checks):
                if sum(candidate_error[v] for v in graph.check_adj[c]) % 2 != syndrome[c]:
                    satisfied = False
                    break
            if satisfied:
                return True, it
        return False, max_iter

    @staticmethod
    def daq_bp_integer(graph: TannerGraph, syndrome: List[int], p: float, max_iter: int = 35) -> Tuple[bool, int]:
        """Our Novel DAQ-BP: 5-bit Integer Quantized, Iteration-Adaptive Damping, Anti-Trapping Momentum."""
        init_llr = 6 # 5-bit channel LLR
        m_c2v = { (c, v): 0 for c in range(graph.num_checks) for v in graph.check_adj[c] }
        m_v2c = { (v, c): init_llr for v in range(graph.num_bits) for c in graph.bit_adj[v] }
        prev_m_v2c = { (v, c): init_llr for v in range(graph.num_bits) for c in graph.bit_adj[v] }
        
        for it in range(1, max_iter + 1):
            alpha_num = 11 if it < 6 else (13 if it < 14 else 14) # Shift divisor = 16
            
            # Check node update (Integer Min-Sum)
            for c in range(graph.num_checks):
                for v in graph.check_adj[c]:
                    sgn = -1 if syndrome[c] == 1 else 1
                    min_abs = 15
                    for v_prime in graph.check_adj[c]:
                        if v_prime == v: continue
                        msg = m_v2c[(v_prime, c)]
                        if msg < 0: sgn *= -1
                        abs_val = abs(msg)
                        if abs_val < min_abs: min_abs = abs_val
                    
                    scaled = (min_abs * alpha_num) >> 4
                    if scaled == 0 and min_abs > 0: scaled = 1
                    m_c2v[(c, v)] = sgn * scaled
            
            # Variable node update with anti-trapping momentum
            candidate_error = [0] * graph.num_bits
            for v in range(graph.num_bits):
                incoming_sum = sum(m_c2v[(c, v)] for c in graph.bit_adj[v])
                total_llr = max(-15, min(15, init_llr + incoming_sum))
                candidate_error[v] = 1 if total_llr < 0 else 0
                
                for c in graph.bit_adj[v]:
                    new_val = init_llr + incoming_sum - m_c2v[(c, v)]
                    prev = prev_m_v2c[(v, c)]
                    if (new_val > 0 and prev < 0) or (new_val < 0 and prev > 0):
                        new_val = new_val >> 1 # Momentum damping
                    
                    prev_m_v2c[(v, c)] = m_v2c[(v, c)]
                    m_v2c[(v, c)] = max(-15, min(15, new_val))
            
            satisfied = True
            for c in range(graph.num_checks):
                if sum(candidate_error[v] for v in graph.check_adj[c]) % 2 != syndrome[c]:
                    satisfied = False
                    break
            if satisfied:
                return True, it
        return False, max_iter

def run_gross_code_sweep():
    H_X, _ = build_gross_code_matrices()
    graph = TannerGraph(H_X) # 72 checks x 144 qubits
    random.seed(2026)
    
    error_rates = [0.01, 0.02, 0.03, 0.04, 0.05]
    trials_per_rate = 200
    
    print(f"=== FULL IBM [[144, 12, 12]] GROSS CODE BENCHMARK ===")
    print(f"Checks: {graph.num_checks}, Physical Qubits: {graph.num_bits}")
    print(f"Trials per rate: {trials_per_rate}\n")
    print(f"{'p_error':<8} | {'Float BP':<12} | {'NMS':<12} | {'Relay-BP (Float)':<18} | {'DAQ-BP (Our Work)':<18} | {'DAQ-BP DSPs'}")
    print("-" * 85)
    
    for p in error_rates:
        succ_bp = succ_nms = succ_relay = succ_daq = 0
        
        for _ in range(trials_per_rate):
            # Physical Pauli-X error vector
            true_error = [1 if random.random() < p else 0 for _ in range(graph.num_bits)]
            if sum(true_error) == 0:
                true_error[random.randint(0, graph.num_bits - 1)] = 1
                
            syndrome = [sum(true_error[v] for v in graph.check_adj[c]) % 2 for c in range(graph.num_checks)]
            
            ok_bp, _ = GrossCodeDecoders.standard_float_bp(graph, syndrome, p)
            if ok_bp: succ_bp += 1
            
            ok_nms, _ = GrossCodeDecoders.classical_nms(graph, syndrome, p)
            if ok_nms: succ_nms += 1
            
            ok_relay, _ = GrossCodeDecoders.relay_bp_float(graph, syndrome, p)
            if ok_relay: succ_relay += 1
            
            ok_daq, _ = GrossCodeDecoders.daq_bp_integer(graph, syndrome, p)
            if ok_daq: succ_daq += 1
            
        rate_bp = (succ_bp / trials_per_rate) * 100
        rate_nms = (succ_nms / trials_per_rate) * 100
        rate_relay = (succ_relay / trials_per_rate) * 100
        rate_daq = (succ_daq / trials_per_rate) * 100
        
        print(f" {p:<7.2f} | {rate_bp:>10.1f}% | {rate_nms:>10.1f}% | {rate_relay:>16.1f}% | {rate_daq:>16.1f}% | 0 DSPs (Integer)")

if __name__ == "__main__":
    run_gross_code_sweep()
