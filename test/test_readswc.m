addpath('../source/');

[A, id, pid, coord, r, subset] = readswc('../data/0-2a.CNG.swc');

fprintf('--- SWC File Summary ---\n');
fprintf('Matrix A:          %d rows × %d cols\n', size(A,1), size(A,2));
fprintf('Number of IDs:     %d\n', numel(id));
fprintf('Number of PIDs:    %d\n', numel(pid));
fprintf('Coordinates:       %d rows × %d cols\n', size(coord,1), size(coord,2));
fprintf('Radii:             %d elements\n', numel(r));
fprintf('Subset labels:     %d entries\n', numel(subset));

% Optionally show unique subset types
uniqueTypes = unique(string(subset));
fprintf('Unique subset types: %s\n', strjoin(uniqueTypes, ', '));