# MATLAB Hodgkin–Huxley Neuron Project

![Demo](output/sbdf2_results/testvideo2.gif)

> Simulate and visualize action potential propagation on real neuron morphologies (SWC), with operator-splitting diffusion–reaction (Hodgkin–Huxley) and publication-ready plots/movies.

---

# MATLAB Hodgkin–Huxley Project

Hodgkin–Huxley (HH) neuron model utilities and solvers in MATLAB, including:
- SWC morphology parsing
- Graph construction & visualization (3D neuron views + Hines sparsity)
- Multiple time-stepping methods (SBDF2; Strang splitting with FE/BE/TR/RK4)
- Trace plotting and movie generation

> This README reflects the repository **as uploaded** (zip contents). Paths and commands match the current `source/`, `data/`, and `output/` layout.

---

## Table of Contents

- [Overview](#overview)
- [Directory Structure](#directory-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Data](#data)
- [Simulation Methods](#simulation-methods)
  - [SBDF2](#sbdf2)
  - [Strang Splitting Variants](#strang-splitting-variants)
- [API Reference](#api-reference)
  - [I/O & Geometry](#io--geometry)
  - [Solvers](#solvers)
  - [Visualization & Postprocessing](#visualization--postprocessing)
  - [Parameters](#parameters)
- [Outputs](#outputs)
- [Reproducibility Tips](#reproducibility-tips)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Overview

This project provides a modular workflow to simulate membrane dynamics on neuron morphologies:

1. **Import morphology** from SWC (`readswc.m`).
2. **Construct graph** (adjacency, neighbors, branch/boundary nodes) and **visualize** (`getgraphstructure.m`, `plotneuron3views.m`).
3. **Assemble stencils** for the cable equation (`stencilmaker.m`).
4. **Integrate in time** using:
   - **SBDF2** (semi-implicit, stiff-friendly) via `sbdf2solve.m`
   - **Strang splitting** with multiple sub-integrators (`strangsolve.m` with `'fe'|'be'|'tr'|'rk4'|...`)
5. **Plot traces** and **export videos** (`plottimeseries.m`, `makematlabmovtraces.m`).

---

## Directory Structure
```
matlabHodgkinHuxleyProject-main
├── data
│   ├── ref1.swc
│   ├── ref2.swc
│   ├── ref3.swc
│   ├── ref4.swc
│   ├── ref5.swc
│   ├── ref6.swc
│   ├── ref7.swc
│   └── ref8.swc
├── output
│   └── sbdf2_results
│       ├── sbdf2neuron.png
│       └── testvideo2.gif
├── source
│   ├── gates.m
│   ├── getgraphstructure.m
│   ├── hh_params.m
│   ├── makematlabmovtraces.m
│   ├── plotneuron3views.m
│   ├── plottimeseries.m
│   ├── readswc.m
│   ├── runsim.m
│   ├── sbdf2solve.m
│   ├── sim_params.m
│   ├── stencilmaker.m
│   ├── strangsolve.m
│   └── timestep.m
├── .gitignore
├── LICENSE
└── README.md
```

- `source/` — MATLAB functions and scripts
- `data/` — example SWC morphologies (`ref*.swc`)
- `output/` — created by scripts at runtime (figures, videos, log files)

---

## Requirements

- **MATLAB** R2019a or newer recommended (uses `readmatrix`, `graph`, etc.)
- Standard toolboxes only (no proprietary toolboxes required)

---

## Installation

From a MATLAB session started at the repository root:

```matlab
addpath('source');        % add project functions to your path
savepath;                 % optional: persist to MATLAB pathdef
```

---

## Quick Start

Run the bundled demo that exercises multiple solvers and exports figures/movies:

```matlab
addpath('source');
run('source/runsim.m');
```

What it does:
- Loads `../data/ref2.swc`
- Plots neuron geometry (3 orthogonal views)
- Runs `sbdf2solve` and several `strangsolve` variants (`'tr'`, `'hn'`, `'be'`, `'fe'`, `'rk4'`)
- Writes results under `output/*`

Minimal example:

```matlab
addpath('source');

% Parse morphology
[A, id, pid, coord, r, subset] = readswc('data/ref2.swc');

% Build graph and visualize
[M, nLst, bLst, brchLst, nNodes, nEdges, meanEdge, maxEdge, minEdge, medEdge] = ...
    getgraphstructure('data/ref2.swc', true, true, false);

% One-step of a solver (example call signature)
dt = 4.0e-5;
sbdf2solve(dt, 1, 105, 'data/ref2.swc', 'output/sbdf2_results', 'SBDF2', true);
```

---

## Data

Example morphologies are included:

- `data/ref1.swc`, `data/ref2.swc`, ..., `data/ref8.swc`

Any valid SWC should work. Header lines (`# ...`) are supported by `readswc.m`.

---

## Simulation Methods

### SBDF2

`sbdf2solve.m` implements **Second-Order Backward Differentiation Formula** with semi-implicit treatment (well-suited to stiff ionic currents / diffusion terms).

**Signature**
```matlab
sbdf2solve(dt, clamp_index, rec_ind, filename, outputFolder, pname, saveall)
```

- `dt` — time step (seconds)
- `clamp_index` — node index for current clamp (stimulus)
- `rec_ind` — receiver index / vector for recording
- `filename` — SWC path (e.g., `'data/ref2.swc'`)
- `outputFolder` — destination for figures/movies
- `pname` — plot/file prefix
- `saveall` — logical flag to write plots/videos

### Strang Splitting Variants

`strangsolve.m` implements **operator splitting** with sub-integrators selectable via the first argument:

- `'fe'` — Forward Euler
- `'be'` — Backward Euler
- `'tr'` — Trapezoid (Crank–Nicolson)
- `'rk4'` — Classical Runge–Kutta 4
- `'hn'` — Heun (explicit trapezoid)

**Signature**
```matlab
strangsolve(method, dt, clamp_index, rec_ind, force, filename, outputFolder, pname, saveall)
```

- `method` — one of `'fe'|'be'|'tr'|'rk4'|'hn'`
- other parameters mirror `sbdf2solve`

---

## API Reference

### I/O & Geometry

- **`readswc.m`**
  ```matlab
  [A, id, pid, coord, r, subset] = readswc(filename)
  ```
  Reads an SWC file (skips `#` headers, spaces/tabs tolerated).  
  Returns numeric matrix `A = [id type x y z r pid]`, integer `id/pid`, coordinates `coord`, radii `r`, and type labels `subset`.

- **`getgraphstructure.m`**
  ```matlab
  [M, nLst, bLst, brchLst, numNodes, numEdges, meanEdge, maxEdge, minEdge, medEdge] = ...
      getgraphstructure(filename, plt, verbose, saveOut)
  ```
  Builds a MATLAB `graph` from SWC, adjacency `M`, neighbor lists, boundary & branch nodes, and edge-length stats.  
  If `plt` or `saveOut` is true, renders 3D neuron and Hines sparsity (`spy`). When `saveOut`, images/logs go to `../output/`.

### Solvers

- **`sbdf2solve.m`** — see [SBDF2](#sbdf2)

- **`strangsolve.m`** — see [Strang Splitting Variants](#strang-splitting-variants)

- **`stencilmaker.m`**
  ```matlab
  A = stencilmaker(n, R, a, C, filename)
  ```
  Assembles system matrices/stencils for the cable equation given node/edge data (`n`, `R`, radii `a`, capacitance `C`, and geometry file).

- **`timestep.m`**
  ```matlab
  result = tr(v, pp, dt, fun, P)   % plus other methods within
  ```
  Time-stepping kernels (e.g., trapezoid, RK, etc.) used by the solvers.

### Visualization & Postprocessing

- **`plotneuron3views.m`** — 3 orthogonal 3D plots of the neuron graph.
  ```matlab
  plotneuron3views(coord, id, pid, outPngPath)
  ```

- **`plottimeseries.m`** — plot time traces for selected nodes/variables.

- **`makematlabmovtraces.m`** — assemble a movie from saved trace images (MP4/GIF depending on your MATLAB/export setup).
  ```matlab
  makematlabmovtraces(dataFolder, geometryFile, outname, idx)
  ```

### Parameters

- **`hh_params.m`**
  ```matlab
  P = hh_params(varargin)
  ```
  Returns a struct of HH parameters (conductances, reversals, capacitance, etc.).

- **`sim_params.m`**
  ```matlab
  S = sim_params(dt, varargin)
  ```
  Returns a struct of simulation controls (time step, duration, clamp, recording indices, seeds, etc.).

- **`gates.m`** — gating kinetics utilities (e.g., `X(v,ss,P)` and sub-functions for α/β rates).

---

## Outputs

Scripts write into `output/` by default, creating subfolders per experiment, e.g.:

- `output/sbdf2_results/*`
- `output/strang0TR_results/*`
- (figures: neuron views, time traces; optional videos via `makematlabmovtraces`)

---

## Reproducibility Tips

- Always call `addpath('source')` in a fresh session.
- Fix the random seed if any stochastic elements are added in future.
- Record the exact `dt`, clamp indices, and SWC used in figure titles or filenames.
- Use `saveOut=true` in `getgraphstructure` to capture logs + figures alongside data.

---

## Troubleshooting

- **`readswc` returns NaNs or zeros**  
  Ensure the SWC numeric lines don’t have unexpected delimiters; the provided parser tolerates leading spaces/tabs and `#` headers.

- **Plots don’t show, but files appear**  
  You likely ran with `saveOut=true` and `plt=false` (off-screen rendering). Set `plt=true` to display.

- **No `.sln` found (Windows build scripts)**  
  The MATLAB project does not require Visual Studio; the PowerShell helper scripts are for UG4 tooling in a separate workflow.

---

## License

See `LICENSE` in the repository root.

---
