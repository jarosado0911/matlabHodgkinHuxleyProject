% MAKE_CONVERGENCE_PLOT_FROM_DATA
% Rebuilds the Figure 2.17A-style SBDF2 convergence plot (voltage at a
% branch point distal from the soma, across refinement levels) from the
% trace_data.mat files already written by run_sbdf2_convergence.m, rather
% than re-running the solver. Used when only some refinement levels of a
% run_sbdf2_convergence.m sweep were allowed to finish (e.g. the finer
% levels were stopped early because of their runtime) but the coarser
% levels' data is still on disk and worth plotting on its own.
%
% Edit `levels` below to match whichever level<k> subfolders actually
% contain a trace_data.mat.

clear; clc;

outBase = '../output/sbdf2_convergence';

levels = struct('folder', {'level0','level1','level2'}, ...
                 'label',  {'0Ref, dt=2.000 [us]','1Ref, dt=1.000 [us]','2Ref, dt=0.500 [us]'});

distalCoord = [-2.490, 785.710, -58.700];  % same physical branch point at every refinement level

figure; hold on;
for k = 1:numel(levels)
    matPath = fullfile(outBase, levels(k).folder, 'trace_data.mat');
    if ~exist(matPath, 'file')
        warning('Skipping missing file: %s', matPath);
        continue;
    end
    S = load(matPath, 't', 'rec_u');
    plot(S.t.*1e3, S.rec_u(1,:).*1e3, 'DisplayName', levels(k).label);
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
