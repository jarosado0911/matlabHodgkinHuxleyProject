function plotstencilheatmaps(LHS, RHS)
    % Choose a common color range from both matrices' nonzeros
    vals = [nonzeros(LHS); nonzeros(RHS)];
    if isempty(vals), vals = 0; end
    clim = [min(vals), max(vals)];
    if clim(1) == clim(2), clim = clim + [-1 1]*eps; end

    % Use tiledlayout if available; otherwise fall back to subplot
    useTiles = exist('tiledlayout','file') == 2;
    if useTiles
        f = figure('Color','w');
        t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        nexttile;  plot_one_matrix(LHS, clim, 'LHS');
        nexttile;  plot_one_matrix(RHS, clim, 'RHS');
        colormap(t, jet);
        cb = colorbar; cb.Layout.Tile = 'east';
    else
        figure('Color','w');
        subplot(1,2,1); plot_one_matrix(LHS, clim, 'LHS');
        subplot(1,2,2); plot_one_matrix(RHS, clim, 'RHS');
        colormap(jet); colorbar;
    end
end

function plot_one_matrix(S, clim, ttl)
    % Sparse-aware: avoid densifying huge matrices
    useSparseMode = issparse(S) && (numel(S) > 2e7);

    if useSparseMode
        [i, j, v] = find(S);
        scatter(j, i, 10, v, 's', 'filled');
        set(gca,'YDir','normal');            % keep same orientation as dense plot
        axis equal tight; xlim([0.5 size(S,2)+0.5]); ylim([0.5 size(S,1)+0.5]);
        caxis(clim); colorbar;
        title(sprintf('%s (sparse view, %d nnz)', ttl, numel(v)), 'Interpreter','none');
    else
        M = full(S);
        h = imagesc(M);                      % create image first
        caxis(clim);
        set(h,'AlphaData', M ~= 0);          % make zeros transparent
        set(gca,'YDir','normal');            % matrix-like orientation
        axis image tight; colorbar;
        title(sprintf('%s (%dx%d)', ttl, size(M,1), size(M,2)), 'Interpreter','none');
    end

    xlabel('column index'); ylabel('row index');
end
