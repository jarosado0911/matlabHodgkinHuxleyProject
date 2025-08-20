# MATLAB Hodgkin–Huxley Neuron Project

![Demo](output/sbdf2_results/testvideo2.gif)

> Simulate and visualize action potential propagation on real neuron morphologies (SWC), with operator-splitting diffusion–reaction (Hodgkin–Huxley) and publication-ready plots/movies.

---

## Table of Contents
- [Overview](#overview)
- [Repository Layout](#repository-layout)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Running a Simulation](#running-a-simulation)
- [Visualizations](#visualizations)
  - [A. 3 views of the neuron](#a-3-views-of-the-neuron)
  - [B. Stencil heatmaps (LHS/RHS)](#b-stencil-heatmaps-lhsrhs)
  - [C. Movie with time-series traces](#c-movie-with-time-series-traces)
  - [D. Plot a single time-series](#d-plot-a-single-time-series)
- [Parameters](#parameters)
  - [Hodgkin–Huxley parameters (`hh_params.m`) †](#hodgkinhuxley-parameters-hh_paramsm-)
  - [Simulation parameters (`sim_params.m`)](#simulation-parameters-sim_paramsm)
- [Data & Outputs](#data--outputs)
- [Convert MP4 → GIF](#convert-mp4--gif)
- [Troubleshooting](#troubleshooting)
- [License](#license)

---

## Overview

This project simulates transmembrane voltage dynamics on a neuron graph extracted from an SWC morphology. The PDE–ODE system is advanced by **Strang operator splitting**: a diffusion linear solve (stencil matrices **LHS/RHS**) alternates with Hodgkin–Huxley **reaction** (gating variable ODEs).

Outputs include:
- Per-node voltage snapshots (`vm_t*.dat`) suitable for movies
- Time vector (`time.mat`)
- Recorded traces at branch nodes (`trace_data.mat`)
- Optional movie (`.mp4`) and demo GIF

---

## Repository Layout

```
data/                         # Sample morphologies (.swc)
output/
  results/
    time.mat                  # time vector t
    trace_data.mat            # rec_u, rec_h, rec_n, rec_m, etc.
  testvideo.mp4               # demo movie
  testvideo2.mp4              # alternate demo movie
  testvideo2.gif              # GIF used in README banner
  *_graph.png, *_sparsity.png # example diagnostics
source/
  getgraphstructure.m
  hh_params.m
  makematlabmov.m
  makematlabmovtraces.m
  neuron_sim.m
  plotneuron3views.m
  plotstencilheatmaps.m
  plottimeseries.m
  print_hh_params.m
  print_sim_params.m
  readswc.m
  sim_params.m
  stencilmaker.m
test/                         # simple tests
README.md
LICENSE
```

---

## Requirements

- **MATLAB** (R2020b+ recommended; base MATLAB only is fine)
- For MP4→GIF (optional):
  - MATLAB-only script provided below, or
  - **FFmpeg** on PATH (Windows: `winget install FFmpeg.FFmpeg`), for best quality

Add the code folder to the MATLAB path once per session:
```matlab
addpath(genpath('source'));
```

---

## Quick Start

Run a short simulation on the included SWC and generate outputs:

```matlab
addpath(genpath('source'));

dt = 1e-5;                              % 10 µs
geom = 'data/0-2a.CNG.swc';             % example morphology
neuron_sim(dt, geom, 1);                % sv=1 saves vm_t*.dat frames

% Optional: movie + traces visualization
makematlabmovtraces('output/results', geom, 'output/testvideo2', 1);
```

This will populate `output/results/` with `time.mat`, `trace_data.mat`, and (if `sv=1`) a `data/` folder of per-frame voltages. The `makematlabmovtraces` call writes `output/testvideo2.mp4` (and you can convert to GIF—see below).

---

## Running a Simulation

```matlab
dt   = 1e-5;                     % time step (seconds)
geom = 'data/0-2a.CNG.swc';      % or data/refinement_*.swc

neuron_sim(dt, geom, 1);         % third arg sv: 1 = save frames, 0 = don't
```
Key behaviors:
- Saves `output/results/time.mat` (`t`) and `output/results/trace_data.mat` (traces)
- If `sv==1`, also saves per-frame voltages as `output/results/data/vm_t<k>.dat`
- The solver prints basic geometry stats and the parameter tables it uses

---

## Visualizations

### A. 3 views of the neuron

```matlab
[A,id,pid,coord,~,~] = readswc('data/0-2a.CNG.swc'); %#ok<ASGLU>
plotneuron3views(coord, id, pid, 'output/neuron_views.png');  % saves PNG
```

### B. Stencil heatmaps (LHS/RHS)

```matlab
% Build matrices and visualize
P = hh_params();
[~,~,~,~,n,~,dx,~,~,~] = getgraphstructure(geom,false,false,false);
[LHS,RHS] = stencilmaker(n, dt, dx, P.R, A(:,6), P.C, geom); % A(:,6) = radii
plotstencilheatmaps(LHS, RHS);  % side-by-side normalized heatmaps
```

### C. Movie with time-series traces

Left: morphology colored by voltage; Right: 4 stacked plots of recorded traces (u/h/n/m) for the selected branch index.

```matlab
makematlabmovtraces('output/results', geom, 'output/testvideo2', 1);  % idx=1
```

> **Tip:** The highlighted large red dot in the left pane shows which **branch node** the right-hand traces correspond to (`idx`th element of the saved `record_index`).

### D. Plot a single time-series

```matlab
S = load('output/results/trace_data.mat');  % loads t, rec_u, rec_h, rec_n, rec_m
plottimeseries(S.t, S.rec_u, 3);            % plot voltage for branch index 3
```

---

## Parameters

### Hodgkin–Huxley parameters (`hh_params.m`) †

```matlab
P = hh_params();                    % defaults
P = hh_params('gna',60e1,'gl',3e1); % override any field(s)
```

Fields (MKS units): `R, C, gk, ek, gna, ena, gl, el`.

† A printed table is shown at runtime via `print_hh_params(P)`.

### Simulation parameters (`sim_params.m`)

```matlab
S = sim_params(1e-5);                              % uses your dt
S = sim_params(1e-5,'endTime',80e-3,'vClamp',60e-3);
```

Fields: `dt, endTime, delay, vStart, vClamp, nT, ni, mi, hi`.  
Printed via `print_sim_params(S)`.

---

## Data & Outputs

- `output/results/time.mat` → `t` vector (seconds)
- `output/results/trace_data.mat` → variables:
  - `t` (T×1), `usoma` (T×1), and branch-node traces:
  - `rec_u` (B×T, volts), `rec_h`, `rec_m`, `rec_n` (B×T, unitless)
  - `record_index` (B×1) — indices of the branch nodes in the morphology
- (optional) `output/results/data/vm_t<k>.dat` → per-node voltage snapshots (volts) for movies

The **branch-node order** in traces matches `record_index`. The movie/trace tools use the same index `idx` to select which branch node to display.

---

## Convert MP4 → GIF

**MATLAB-only** conversion (no FFmpeg needed):

```matlab
function mp4_to_gif(mp4Path, gifPath, fps, maxDim)
if nargin < 3, fps = 12; end
if nargin < 4, maxDim = 640; end
v = VideoReader(mp4Path); dt = 1/fps; nextT = 0; first = true;
while hasFrame(v)
    f = readFrame(v);
    if v.CurrentTime+1e-6 < nextT, continue; end, nextT = nextT + dt;
    scale = min(1, maxDim / max(size(f,1), size(f,2))); if scale<1, f = imresize(f, scale); end
    if first, [A,map] = rgb2ind(f,256,'nodither'); imwrite(A,map,gifPath,'gif','LoopCount',Inf,'DelayTime',dt); first=false;
    else, A = rgb2ind(f,map,'nodither'); imwrite(A,map,gifPath,'gif','WriteMode','append','DelayTime',dt); end
end
end
```
Usage:
```matlab
mp4_to_gif('output/testvideo2.mp4','output/testvideo2.gif',12,640);
```

**FFmpeg** (best quality, if installed):
```bash
ffmpeg -y -i output/testvideo2.mp4 -vf "fps=12,scale=640:-1:flags=lanczos,palettegen" /tmp/pal.png
ffmpeg -y -i output/testvideo2.mp4 -i /tmp/pal.png -lavfi "fps=12,scale=640:-1:flags=lanczos,paletteuse=dither=bayer:bayer_scale=5" output/testvideo2.gif
```

---

## Troubleshooting

- **Blank/white figure on save:** set figure `'InvertHardcopy','off'` (already done).
- **Minor “subgrid” lines visible:** in axes, set `XMinorGrid/YMinorGrid/ZMinorGrid` to `'off'` (done in `style_graypanel`).
- **Colorbar placement inside axes:** we manually position with `'Location','manual'` and `InnerPosition`.
- **Edges color don’t match nodes:** ensure you use `makematlabmovtraces` (interpolated edge colors from endpoints).
- **Git ignoring generated frames:** add to `.gitignore`
  ```
  output/results/data/*
  !output/results/data/.gitkeep
  ```

---

## License

This project is provided under the terms in `LICENSE` (MIT-style). Please cite appropriately if you use this code in your work.
