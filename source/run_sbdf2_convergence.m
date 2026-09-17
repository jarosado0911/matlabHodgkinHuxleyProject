% RUN_SBDF2_CONVERGENCE
% Reproduces Figure 2.17A from the dissertation: the SBDF2 spatial/temporal
% convergence of the membrane voltage trace at a branch point distal from
% the soma, as Delta x, Delta t -> 0. Unlike the dissertation figure, this
% does NOT include the NEURON reference curve - MatLab/SBDF2 only, using
% the Rank-1 (Sherman-Morrison) voltage-clamp update already implemented
% in sbdf2solve.m/DircheletRank1UpdateSolve.
%
% Refinement levels: ref3.swc..ref7.swc were verified to be the same
% morphology at successive edge-splitting refinements, with mean edge
% length 4, 2, 1, 0.5, 0.25 [um] respectively - i.e. the dissertation's
% 0th..4th refinement - each paired with the matching halved time step
% (2, 1, 0.5, 0.25, 0.125 [us]), exactly as in Figure 2.17.
%
% The recording site is chosen automatically as the branch point (degree
% > 2 in the neuron graph) that is farthest (Euclidean distance) from the
% soma - the same physical branch point at every refinement level, even
% though its node index changes as edges are subdivided.

clear; clc;

files  = {'ref3.swc','ref4.swc','ref5.swc','ref6.swc','ref7.swc'};
dts    = [2.0e-6, 1.0e-6, 0.5e-6, 0.25e-6, 0.125e-6];
labels = {'0Ref, dt=2.000 [us]','1Ref, dt=1.000 [us]','2Ref, dt=0.500 [us]', ...
          '3Ref, dt=0.250 [us]','4Ref, dt=0.125 [us]'};

dataDir = '../data';
outBase = '../output/sbdf2_convergence';
if ~exist(outBase,'dir'), mkdir(outBase); end

clamp_index = 1;  % soma is always node id 1 in these regularized/refined SWCs

figure; hold on;

distalCoord = [];
for k = 1:numel(files)
    swcfile = fullfile(dataDir, files{k});

    [rec_ind, distalCoord, distalDist] = find_distal_branch(swcfile);

    fprintf('%s: recording at branch node %d, coord = (%.3f, %.3f, %.3f), distance from soma = %.3f [um]\n', ...
        files{k}, rec_ind, distalCoord(1), distalCoord(2), distalCoord(3), distalDist);

    levelOut = fullfile(outBase, sprintf('level%d', k-1));

    % saveall = false: skip the per-timestep snapshot dump (not needed for
    % this trace-only convergence plot) - the recorded trace itself is
    % still written to <levelOut>/trace_data.mat by sbdf2solve.m.
    sbdf2solve(dts(k), clamp_index, rec_ind, swcfile, levelOut, labels{k}, false);
end

ylim([-10, 50]);
title(sprintf('Voltage at distal branch point = (%.3f, %.3f, %.3f)', ...
    distalCoord(1), distalCoord(2), distalCoord(3)));
xlabel('time [ms]');
ylabel('voltage [mV]');
legend('show', 'Location', 'best');

pngPath = fullfile(outBase, 'sbdf2_convergence_branchA.png');
figPath = fullfile(outBase, 'sbdf2_convergence_branchA.fig');
exportgraphics(gcf, pngPath, 'Resolution', 200);
savefig(gcf, figPath);

fprintf('Saved figure to %s\n', pngPath);

% ------------------------------------------------------------------------
function [distalId, distalCoord, distalDist] = find_distal_branch(filename)
% Finds the branch node (degree > 2 in the undirected neuron graph) that
% is farthest (Euclidean distance) from the soma (node id 1).
    [~,~,~,brchLst,~,~,~,~,~,~] = getgraphstructure(filename, false, false, false);
    [~,~,~,coord,~,~] = readswc(filename);

    somaCoord = coord(1,:);
    nB = numel(brchLst);
    dists = zeros(nB,1);
    for i = 1:nB
        dists(i) = norm(coord(brchLst{i},:) - somaCoord);
    end
    [distalDist, idx] = max(dists);
    distalId = brchLst{idx};
    distalCoord = coord(distalId,:);
end
