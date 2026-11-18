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

## Experimental Results (Preliminary)

Syndrome extraction memory experiment on the IBM [[144, 12, 12]] Gross code across 3 rounds under depolarizing noise (1,000 shots per point):

| Physical Noise $p$ | Float-BP LER | NMS LER ($\alpha=0.75$) | Relay-BP LER ($\gamma=-0.25$) | **DAQ-BP (5-bit Integer)** | Relative LER Gain vs Float-BP |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 0.001 | 2.30% | 2.10% | 3.10% | **1.40%** | **-39.1%** |
| 0.002 | 4.50% | 4.50% | 7.50% | **3.70%** | **-17.8%** |
| 0.004 | 8.60% | 6.30% | 9.90% | **5.30%** | **-38.4%** |
| 0.006 | 12.50% | 11.00% | 11.80% | **8.90%** | **-28.8%** |
| 0.008 | 15.40% | 13.20% | 18.80% | **9.90%** | **-35.7%** |
| 0.010 | 18.90% | 18.30% | 20.10% | **14.40%** | **-23.8%** |

*Note: Results reflect accumulated syndrome decoding over spatial checks. Relative LER reductions range from 17.8% to 39.1% across tested noise rates. Future work includes evaluating against full 3D space-time detector error models (DEM) and BP-OSD baselines with $\ge 10^5$ shots.*

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

