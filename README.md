# DAQ-BP: Degeneracy-Aware Quantized Belief Propagation Decoder

> **Status: Research Prototype / Preliminary Work**  
> This repository contains the reference implementation, experimental benchmarks, and draft manuscript for DAQ-BP, an integer-quantized belief propagation decoder for quantum LDPC codes.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

Hardware-efficient, integer-only belief propagation decoder designed for quantum low-density parity-check (qLDPC) codes. Evaluated on the IBM [[144, 12, 12]] bivariate bicycle Gross code using Google's Stim simulator.

## Key Features

- **5-bit Integer Quantization**: Complete elimination of 64-bit floating-point transcendental functions and DSP multipliers.
- **Iteration-Adaptive Damping**: Dynamic bit-shift attenuation ($\gg \lfloor k/10 \rfloor$) to mitigate short-cycle correlation buildup.
- **Anti-Trapping Momentum**: Sign-inversion tracking that halves message momentum upon oscillation to escape degenerate trapping sets.
- **Pipelined RTL Architecture**: Synthesizable SystemVerilog models for sub-60 ns syndrome resolution on modern FPGAs (analytical timing model).

## Experimental Results (Code-Capacity Model)

Evaluated on the IBM [[144, 12, 12]] Gross code under independent data-qubit Pauli-X noise across 2,000 Monte Carlo shots per point:

### 1. Logical Error Rate (LER)

| Physical Error $p$ | Float-BP LER | NMS LER ($\alpha=0.75$) | Relay-BP LER | **DAQ-BP LER (5-bit Integer)** |
| :---: | :---: | :---: | :---: | :---: |
| 0.005 | 0.35% | 0.15% | 0.40% | **0.30%** |
| 0.010 | 1.25% | 0.90% | 1.55% | **1.25%** |
| 0.020 | 4.75% | 3.00% | 6.95% | **4.70%** |
| 0.030 | 10.35% | 8.80% | 13.95% | **10.70%** |
| 0.040 | 13.10% | 12.55% | 17.50% | **15.15%** |
| 0.050 | 19.50% | 17.85% | 23.65% | **22.40%** |

### 2. Convergence Dynamics & Latency Breakdown (DAQ-BP on Non-Zero Syndromes)

| Physical Error $p$ | Non-Zero Conv. % | Mean Iters (Converged) | Mean Latency (Converged) | Timeout % (35 Iters) | Worst-Case Tail Latency |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 0.005 | **96.9%** | **1.04** | **10.4 ns** | 3.1% | 350 ns |
| 0.010 | **90.8%** | **1.06** | **10.6 ns** | 9.2% | 350 ns |
| 0.020 | **70.2%** | **1.15** | **11.5 ns** | 29.8% | 350 ns |
| 0.030 | 43.9% | 1.16 | 11.6 ns | 56.1% | 350 ns |
| 0.040 | 27.0% | 1.18 | 11.8 ns | 73.0% | 350 ns |
| 0.050 | 13.2% | 1.10 | 11.0 ns | 86.8% | 350 ns |

*Note: Latency is based on an estimated 4-stage pipeline running at 400 MHz (2.5 ns clock, 10 ns per full iteration). At physical error rates $p \le 1\%$, DAQ-BP resolves over 90% of non-zero syndromes in ~10.5 ns, with worst-case tail latency bounded at 350 ns (35 iterations).*

## Directory Layout

```text
├── src/                # Core Python package (Gross code constructor & decoders)
├── rtl/                # Synthesizable SystemVerilog implementation
├── benchmarks/         # Multi-process simulation scripts
└── paper/              # LaTeX manuscript and compiled PDF
```

## Quickstart

```bash
git clone https://github.com/Deeven-Seru/daq-bp.git
cd daq-bp
pip install -e .
python benchmarks/run_stim_circuit_level.py
```

## Hardware Model Notes

The hardware latency budget is based on a 4-stage pipelined check-and-variable node processing unit running at $400\,\text{MHz}$ ($2.5\,\text{ns}$ clock period, $10\,\text{ns}$ per complete BP iteration). At an average convergence depth of $4.4$ to $5.5$ iterations, the analytical latency is $44\,\text{ns} - 55\,\text{ns}$. Physical place-and-route reports in AMD Vivado are in preparation.

## License

MIT License. See `LICENSE` for details.

