%% Stencil Heatmap Demo
% Visualize the left- and right-hand stencil matrices produced by your
% discretization using a normalized ([0,1]) jet colormap.
%
% This script:
%   1) Loads a neuron geometry from an SWC file.
%   2) Converts geometric quantities to SI units (meters).
%   3) Builds Hodgkin–Huxley parameter struct (P).
%   4) Retrieves graph information (node count n, average edge length dx).
%   5) Constructs the diffusion/reaction stencil matrices (LHS, RHS).
%   6) Plots both matrices side-by-side as heatmaps with a shared colorbar.
%
% Requirements (on MATLAB path):
%   - readswc.m             : SWC reader returning id/pid/coords/radius/subsets
%   - getgraphstructure.m   : Graph/topology extraction from SWC (returns n, dx, ...)
%   - stencilmaker.m        : Builds LHS and RHS stencil matrices
%   - hh_params.m           : Returns Hodgkin–Huxley parameters (struct P)
%   - plotstencilheatmaps.m : Heatmap visualizer that normalizes to [0,1]
%
% Output:
%   - A figure with two heatmaps (LHS, RHS) using the same color range [0,1];
%     structural zeros appear as blue.
%
% Notes:
%   - Ensure consistent units: radii and spatial step (dx) should be in meters,
%     voltages in volts, time in seconds.
%   - If your getgraphstructure returns dx in micrometers (µm), convert to meters.

%% --- 1) Load geometry (SWC) ---
filename = '../data/0-2a.CNG.swc';
[A, id, pid, coord, r, subset] = readswc(filename); %#ok<ASGLU> % A/coord/subset may be unused downstream

%% --- 2) Units: convert radii to meters (SWC radii typically in micrometers) ---
r = r * 1e-6;   % [µm] -> [m]

%% --- 3) Time step & HH parameters ---
dt = 1e-5;      % [s] time step used by the stencil
P  = hh_params();  % returns struct with fields: R, C, gk, ek, gna, ena, gl, el (SI units)

%% --- 4) Graph information (node count n, average edge length dx) ---
% The exact outputs of getgraphstructure may vary with your implementation.
% Here we only retain n (number of nodes) and dx (average edge length).
[~,~,~,~, n, ~, dx, ~, ~, ~] = getgraphstructure(filename, false, false, false);

% IMPORTANT: If dx is in micrometers (µm), uncomment the next line:
dx = dx * 1e-6;  % [µm] -> [m]

%% --- 5) Build stencil matrices ---
% stencilmaker signature (as used here):
%   [LHS, RHS] = stencilmaker(n, dt, dx, P.R, r, P.C, filename)
% where:
%   n    = number of nodes
%   dt   = time step [s]
%   dx   = representative spatial step [m] (average; per-edge dx usually handled internally)
%   P.R  = axial resistivity [ohm·m]
%   r    = per-node radius vector [m]
%   P.C  = membrane capacitance [F/m^2]
%   filename can be used internally by your builder (e.g., boundary sets)
[LHS, RHS] = stencilmaker(n, dt, dx, P.R, r, P.C, filename);

%% --- 6) Visualize stencils as normalized heatmaps ---
% plotstencilheatmaps(LHS, RHS) will:
%   - gather nonzero entries from both matrices,
%   - linearly scale them into [0,1] (structural zeros stay 0),
%   - show LHS and RHS side-by-side with a shared jet colorbar (0→1).
plotstencilheatmaps(LHS, RHS);
