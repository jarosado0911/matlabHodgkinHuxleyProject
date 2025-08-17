function plotstencilheatmaps(LHS, RHS)
%PLOTSTENCILHEATMAPS  Visualize two stencil matrices as [0,1]-scaled heatmaps.
%   PLOTSTENCILHEATMAPS(LHS, RHS) plots the left-hand-side (LHS) and
%   right-hand-side (RHS) matrices side-by-side using a shared color scale
%   in the range [0, 1]. Structural zeros remain exactly 0 and therefore
%   render as the blue end of the 'jet' colormap. Both axes are forced to
%   the same on-screen size with a single shared colorbar on the right.
%
%   Inputs
%   ------
%   LHS : numeric matrix (dense or sparse)
%       Stencil matrix (e.g., from diffusion solve).
%   RHS : numeric matrix (dense or sparse)
%       Stencil matrix to compare against LHS.
%
%   Behavior
%   --------
%   * Both matrices are linearly normalized to [0,1] using the global min
%     and max of their *nonzero* entries. Structural zeros remain 0.
%   * A single colorbar (0→1) is shared by both plots.
%   * Axes are positioned explicitly so both images have identical size.
%
%   Example
%   -------
%       [LHS,RHS] = stencilmaker(...);
%       plotstencilheatmaps(LHS, RHS);
%
%   Notes
%   -----
%   * For very large sparse matrices, consider creating a sparse-aware
%     plotting path (e.g., scatter of nonzeros) if rendering is slow.

    % ---- compute global scaling using nonzero values from BOTH matrices ----
    nz = [nonzeros(LHS); nonzeros(RHS)];
    if isempty(nz), nz = 0; end
    vmin = min(nz);  % minimum over nonzero entries
    vmax = max(nz);  % maximum over nonzero entries
    if vmin == vmax
        % Avoid degenerate scale: expand slightly so imagesc has range
        vmin = vmin - 1; 
        vmax = vmax + 1;
    end

    % Scale each matrix to [0,1], keeping structural zeros at 0
    LHSs = scale01_keep_zeros(LHS, vmin, vmax);
    RHSs = scale01_keep_zeros(RHS, vmin, vmax);

    % ---- figure + two axes with identical size/position ----
    figure('Color','w');
    ax1 = axes('Units','normalized','Position',[0.08 0.11 0.38 0.80]); % left panel
    ax2 = axes('Units','normalized','Position',[0.54 0.11 0.38 0.80]); % right panel

    % Plot with shared CLim [0,1]; zeros (structural) appear as blue in 'jet'
    imagesc(ax1, LHSs, [0 1]); 
    set(ax1,'YDir','normal'); axis(ax1,'image'); title(ax1,'LHS');

    imagesc(ax2, RHSs, [0 1]); 
    set(ax2,'YDir','normal'); axis(ax2,'image'); title(ax2,'RHS');

    % Shared colormap and a single colorbar on the far right
    colormap(jet);
    cb = colorbar(ax2);
    set(cb, 'Units','normalized', 'Position',[0.94 0.11 0.02 0.80]);
    cb.Label.String = 'Normalized value (0 \rightarrow 1)';

    % Tidy labels and identical formatting
    xlabel(ax1,'column'); ylabel(ax1,'row');
    xlabel(ax2,'column'); ylabel(ax2,'row');
end

function Ms = scale01_keep_zeros(M, vmin, vmax)
%SCALE01_KEEP_ZEROS  Normalize nonzeros of M to [0,1], keep zeros at 0.
%   Ms = SCALE01_KEEP_ZEROS(M, vmin, vmax) linearly maps the nonzero
%   entries of M from [vmin, vmax] to [0,1]. Structural zeros remain 0.
%
%   Inputs
%     M     : numeric matrix (dense or sparse)
%     vmin  : lower bound used for scaling (scalar)
%     vmax  : upper bound used for scaling (scalar)
%
%   Output
%     Ms    : double matrix the same size as M, with values in [0,1]
%             and structural zeros preserved as 0.

    % Initialize output; ensure a floating type for imagesc
    Ms = zeros(size(M), 'like', full(0));

    if vmin == vmax
        % Degenerate case already guarded above; nothing to scale
        return; 
    end

    % Scale only nonzero entries; preserve structural zeros at 0 (blue)
    idx = (M ~= 0);
    Ms(idx) = (full(M(idx)) - vmin) / (vmax - vmin);

    % Clamp numerically to [0,1] in case of round-off
    Ms(Ms < 0) = 0; 
    Ms(Ms > 1) = 1;
end