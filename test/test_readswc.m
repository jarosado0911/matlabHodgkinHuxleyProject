% Add the '../source/' directory to the MATLAB search path.
% This ensures MATLAB can find custom functions like readswc.m.
addpath('../source/');  

% Call readswc to parse the SWC file containing neuron morphology.
% Outputs:
%   A      : full numeric matrix [id, type, x, y, z, r, pid]
%   id     : node IDs
%   pid    : parent IDs (connectivity)
%   coord  : [x,y,z] coordinates of each node
%   r      : radii of each node
%   subset : cell array of node type labels (e.g., 'soma', 'axon', 'dendrite')
[A, id, pid, coord, r, subset] = readswc('../data/0-2a.CNG.swc');

% Print a summary of the SWC file contents.
fprintf('--- SWC File Summary ---\n');
fprintf('Matrix A:          %d rows × %d cols\n', size(A,1), size(A,2));
fprintf('Number of IDs:     %d\n', numel(id));
fprintf('Number of PIDs:    %d\n', numel(pid));
fprintf('Coordinates:       %d rows × %d cols\n', size(coord,1), size(coord,2));
fprintf('Radii:             %d elements\n', numel(r));
fprintf('Subset labels:     %d entries\n', numel(subset));

% Optionally, display the unique node types present in the neuron.
% Converts 'subset' to string array, finds unique values, then joins into a comma-separated list.
uniqueTypes = unique(string(subset));
fprintf('Unique subset types: %s\n', strjoin(uniqueTypes, ', '));