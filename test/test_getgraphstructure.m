% Add the '../source/' directory to the MATLAB search path.
% This allows MATLAB to find your custom functions (like readswc.m, getgraphstructure.m, etc.)
% even if you call them from outside that folder.
addpath('../source/');  

% Specify the SWC neuron morphology file to analyze.
% SWC files contain neuronal structures as node coordinates and connectivity.
filename = '../data/0-2a.CNG.swc';

% Set runtime options:
plt      = true;   % If true, display plots on screen (3D neuron geometry + sparsity pattern).
verbose  = true;   % If true, print detailed information about nodes and edges to the console/log.
saveout  = true;   % If true, save all verbose output to ../output/neuron.log and plots as PNG images.

% Call getgraphstructure to construct the graph from the SWC file.
% The function returns:
%   M         : adjacency matrix (sparse)
%   nLst      : neighbor list for each node (cell array)
%   bLst      : list of boundary nodes (endpoints)
%   brchLst   : list of branching nodes (bifurcations)
%   numNodes  : total number of nodes
%   numEdges  : total number of edges
%   meanEdge  : average edge length (microns)
%   maxEdge   : maximum edge length
%   minEdge   : minimum edge length
%   medEdge   : median edge length
%
% Behavior controlled by the flags:
%   - If plt == true, plots are shown on screen.
%   - If saveout == true, log and plot images are saved into ../output/.
%   - If verbose == true, detailed text output is printed (and logged if saveout == true).
[M,nLst,bLst,brchLst,numNodes,numEdges,meanEdge,maxEdge,minEdge,medEdge] = ...
    getgraphstructure(filename, plt, verbose, saveout);
