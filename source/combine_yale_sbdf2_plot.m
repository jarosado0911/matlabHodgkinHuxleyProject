% COMBINE_YALE_SBDF2_PLOT
% Overlays the NEURON reference trace (yale_neuron_trace.fig, produced by
% runsim_yale_neuron.m) on top of the SBDF2 refinement-level traces
% (sbdf2_convergence_branchA.fig, produced by run_sbdf2_convergence.m /
% make_convergence_plot_from_data.m) into a single combined figure. Both
% figures record at the same physical branch point, so this is a direct
% NEURON-vs-SBDF2 convergence comparison, not just a visual overlay.
%
% Edit yaleFig/sbdf2Fig below if your files live elsewhere.

clear; clc;

yaleFig  = '../output/yale_neuron_results/yale_neuron_trace.fig';
sbdf2Fig = '../output/sbdf2_convergence/sbdf2_convergence_branchA.fig';

outfolder = '../output/combined_results';
if ~exist(outfolder,'dir'), mkdir(outfolder); end

% ---- pull the line data out of each source .fig -------------------------
hYale  = openfig(yaleFig,  'invisible');
hSbdf2 = openfig(sbdf2Fig, 'invisible');

axYale  = findobj(hYale,  'Type', 'axes');
axSbdf2 = findobj(hSbdf2, 'Type', 'axes');

lineYale  = findobj(axYale,  'Type', 'line');
linesSbdf2 = findobj(axSbdf2, 'Type', 'line');

titleText = get(get(axSbdf2, 'Title'), 'String');

% ---- build the combined figure ------------------------------------------
figure; ax = gca; hold(ax, 'on');

copyobj(linesSbdf2, ax);
neuronLine = copyobj(lineYale, ax);
set(neuronLine, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);

close(hYale); close(hSbdf2);

ylim(ax, [-10, 50]);
title(ax, titleText);
xlabel(ax, 'time [ms]');
ylabel(ax, 'voltage [mV]');
legend(ax, 'show', 'Location', 'best');

pngPath = fullfile(outfolder, 'combined_yale_sbdf2.png');
figPath = fullfile(outfolder, 'combined_yale_sbdf2.fig');
exportgraphics(gcf, pngPath, 'Resolution', 200);
savefig(gcf, figPath);

fprintf('Saved figure to %s\n', pngPath);
