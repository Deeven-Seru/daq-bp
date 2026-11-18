# DAQ-BP: Degeneracy-Aware Quantized Belief Propagation Decoder

> **Status: Research Prototype / Preliminary Work**  
> This repository contains the reference implementation, experimental benchmarks, and draft manuscript for DAQ-BP, an integer-quantized belief propagation decoder for quantum LDPC codes.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24333038.svg)](https://doi.org/10.5281/zenodo.24333038)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

Hardware-efficient, integer-only belief propagation decoder designed for quantum low-density parity-check (qLDPC) codes. Evaluated on the IBM [[144, 12, 12]] bivariate bicycle Gross code using Google's Stim simulator.

## Key Features

- **5-bit Integer Quantization**: Complete elimination of 64-bit floating-point transcendental functions and DSP multipliers.
- **Iteration-Adaptive Damping**: Dynamic bit-shift attenuation ($\gg \lfloor k/10 \rfloor$) to mitigate short-cycle correlation buildup.
- **Anti-Trapping Momentum**: Sign-inversion tracking that halves message momentum upon oscillation to escape degenerate trapping sets.
- **Pipelined RTL Architecture**: Synthesizable SystemVerilog models for sub-60 ns syndrome resolution on modern FPGAs (analytical timing model).

## Key Empirical Findings (With Standard BP-OSD Reference)

Evaluated on the IBM [[144, 12, 12]] Gross code against standard BP-OSD (OSD-CS order 10 from the `ldpc` library) across 1,000 Monte Carlo shots per point:

### 1. Logical Error Rate Comparison

| Physical Error $p$ | Float-BP LER | NMS LER | Relay-BP LER | **DAQ-BP (5-bit)** | **BP-OSD (Float Reference)** | **DAQ-BP+OSD (Hybrid)** |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0.005 | 0.20% | 0.10% | 0.20% | **0.10%** | 0.00% | **0.00%** |
| 0.010 | 1.40% | 0.80% | 1.30% | **1.20%** | 0.00% | **0.00%** |
| 0.020 | 4.40% | 2.60% | 6.30% | **4.30%** | 0.20% | **0.20%** |
| 0.030 | 10.60% | 9.20% | 14.60% | **11.10%** | 0.40% | **0.40%** |
| 0.040 | 13.40% | 12.30% | 16.30% | **14.20%** | 1.00% | **1.00%** |
| 0.050 | 19.30% | 18.30% | 21.70% | **21.40%** | 2.30% | **2.30%** |

### 2. Elimination of $O(N^3)$ Matrix Inversions

| Physical Error $p$ | DAQ-BP Non-Zero Conv. % | Mean Converged Iters | OSD Fallback Rate | **Matrix Inversions Avoided** |
| :---: | :---: | :---: | :---: | :---: |
| 0.005 | **97.7%** | **1.03** | **1.2%** | **98.8%** |
| 0.010 | **90.7%** | **1.07** | **6.9%** | **93.1%** |
| 0.020 | **71.6%** | **1.16** | 27.1% | **72.9%** |
| 0.030 | 43.9% | 1.15 | 55.4% | 44.6% |
| 0.040 | 26.0% | 1.12 | 73.9% | 26.1% |
| 0.050 | 13.3% | 1.10 | 86.6% | 13.4% |

*Key Systems Result: Standard BP-OSD invokes $O(N^3)$ Gaussian elimination on 100% of syndrome extraction cycles. By deploying 5-bit integer DAQ-BP as an ultra-fast hardware front-end ($10.5\,\text{ns}$ latency, 0 DSP multipliers), **$93.1\%$ to $98.8\%$ of all matrix inversions are completely eliminated** while maintaining identical logical error suppression.*

## Interactive Validation Notebook

An interactive, executable demonstration and visual verification suite is available under `notebooks/`:
- **Jupyter Notebook:** [`notebooks/daq_bp_validation.ipynb`](notebooks/daq_bp_validation.ipynb)
- **Rendered PDF Report:** [`notebooks/daq_bp_interactive_validation.pdf`](notebooks/daq_bp_interactive_validation.pdf)

The notebook provides:
1. Matrix sparsity visualization of the IBM $[[144, 12, 12]]$ check matrices $H_X$ and $H_Z$.
2. Step-by-step 5-bit lookup table plots comparing quantized vs. continuous hyperbolic functions.
3. Live multi-qubit Pauli-X error injection walkthrough.
4. Comparative Monte Carlo plots illustrating LER pseudo-thresholds and the elimination of $>93\%$ of $O(N^3)$ matrix inversions.

## Directory Layout

```text
├── src/                # Core Python package (Gross code constructor & decoders)
├── rtl/                # Synthesizable SystemVerilog implementation
├── benchmarks/         # Multi-process simulation scripts (including BP-OSD)
├── notebooks/          # Executable Jupyter notebook & rendered PDF report
└── paper/              # LaTeX manuscript and compiled PDF
```

## Quickstart

```bash
git clone https://github.com/Deeven-Seru/daq-bp.git
cd daq-bp
pip install -e .
python benchmarks/benchmark_code_capacity.py
```

## Hardware Model Notes

The hardware latency budget is based on a 4-stage pipelined check-and-variable node processing unit running at $400\,\text{MHz}$ ($2.5\,\text{ns}$ clock period, $10\,\text{ns}$ per complete BP iteration). At operational error rates ($p \le 1\%$), average converged depth is $1.03 - 1.07$ iterations, yielding an analytical latency of $\sim 10.5\,\text{ns}$ and a strict tail latency bound of $350\,\text{ns}$ (35 iterations). Physical place-and-route reports in AMD Vivado are in preparation.

## Citation

If you use this work or benchmark data, please cite:

```bibtex
@misc{seru2026daqbp,
  author       = {Deeven Seru},
  title        = {{DAQ-BP: Degeneracy-Aware Quantized Belief Propagation for Real-Time Decoding of Quantum LDPC Codes}},
  month        = nov,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.24333038},
  url          = {https://doi.org/10.5281/zenodo.24333038}
}
```

## License

MIT License. See `LICENSE` for details.

