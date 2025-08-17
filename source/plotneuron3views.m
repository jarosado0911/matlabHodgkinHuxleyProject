function plotneuron3views(coord, id, pid, outFile)
% plotneuron3views_distcolor(coord,id,pid[,outFile])
% Left: XY (big). Top-right: XZ. Bottom-right: YZ.
% Lines & points colored by distance from origin (jet).
% Panels use a gray background with light gridlines (screenshot style).
% If outFile is provided (e.g., 'figs/neuron.png' or '.pdf'), the figure is saved.
%
% Example:
%   [A,id,pid,coord,~,~] = readSWC('cell.swc'); %#ok<NASGU>
%   plotneuron3views_distcolor(coord,id,pid,'figs/neuron_views.png');

    % ---- color by distance from origin ----
    d = sqrt(sum(coord.^2, 2));
    dmax = max(d); if dmax == 0, dmax = 1; end

    % ---- figure ----
    f = figure('Color',[0 0 0],'InvertHardcopy','off'); % keep dark bg when exporting

    % ---- layout knobs (normalized units) ----
    L        = 0.08;   % outer left margin
    R        = 0.14;   % outer right margin (increase for more right-side space)
    B        = 0.10;   % bottom margin
    T        = 0.08;   % top margin (reserved)
    cbW      = 0.025;  % colorbar width
    gapCB    = 0.045;  % gap between colorbar and LEFT axes (increase to push CB away)
    gapLR    = 0.005;  % gap between LEFT column and RIGHT column (decrease to pull closer)
    leftFrac = 0.56;   % fraction of (available width) given to left panel

    % Total width available after margins + colorbar
    availW = 1 - L - R - cbW - gapCB;
    leftW  = leftFrac * availW;
    rightW = availW - leftW - gapLR;

    % Positions
    xLeft  = L + cbW + gapCB;
    xRight = xLeft + leftW + gapLR;

    axLpos  = [xLeft,  B,     leftW, 0.82]; % LEFT (XY)
    axTRpos = [xRight, 0.56,  rightW, 0.36]; % XZ (top-right)
    axBRpos = [xRight, 0.12,  rightW, 0.36]; % YZ (bottom-right)

    % Build axes (innerposition constraint keeps drawable areas stable)
    axL  = axes('Parent',f,'Units','normalized','Position',axLpos, ...
                'PositionConstraint','innerposition');
    axTR = axes('Parent',f,'Units','normalized','Position',axTRpos, ...
                'PositionConstraint','innerposition');
    axBR = axes('Parent',f,'Units','normalized','Position',axBRpos, ...
                'PositionConstraint','innerposition');

    % ---- draw neuron (distance-colored) ----
    draw_swc_colored(axL,  coord, id, pid, d, 1.2); view(axL,  0, 90); % XY
    draw_swc_colored(axTR, coord, id, pid, d, 1.2); view(axTR, 0,  0); % XZ
    draw_swc_colored(axBR, coord, id, pid, d, 1.2); view(axBR, 90, 0); % YZ

    % ---- enforce equal aspect & identical limits across all panels ----
    set_equal_cube_axes([axL, axTR, axBR], coord, 0.05); % 5% padding
    for ax = [axL, axTR, axBR]
        daspect(ax,[1 1 1]); pbaspect(ax,[1 1 1]); axis(ax,'vis3d');
    end

    % ---- panel styling (gray with gridlines) ----
    apply_graygrid(axL);
    apply_graygrid(axTR);
    apply_graygrid(axBR);

    % Titles (put left title above the axes) & labels
    title(axL,'XY (top)','Color','w','FontWeight','bold', ...
          'Units','normalized','Position',[0.5, 1.03, 0]);
    title(axTR,'XZ (side)','Color','w','FontWeight','bold');
    title(axBR,'YZ (side)','Color','w','FontWeight','bold');

    xlabel(axL,'x','Color','w');  ylabel(axL,'y','Color','w');
    xlabel(axTR,'x','Color','w'); ylabel(axTR,'z','Color','w');
    xlabel(axBR,'y','Color','w'); ylabel(axBR,'z','Color','w');

    % ---- colormap & shared range ----
    colormap(f, jet);
    for ax = [axL, axTR, axBR], caxis(ax,[0 dmax]); end

    % ---- ensure right two panels have identical inner widths ----
    match_right_axes_width(axTR, axBR);

    % ---- colorbar on the LEFT of the left axes ----
    cb = colorbar(axL,'Location','westoutside');
    cb.Color = 'w';
    cb.Label.String = 'Distance from origin';
    cb.Label.Color  = 'w';
    cb.Position = [L, axLpos(2), cbW, axLpos(4)]; % [x y w h]; sits to the left

    % ---- optional save ----
    if nargin >= 4 && ~isempty(outFile)
        [outDir,~,~] = fileparts(outFile);
        if ~isempty(outDir) && ~exist(outDir,'dir'), mkdir(outDir); end
        try
            % Robust export with a bit of padding so nothing gets clipped
            exportgraphics(gcf, outFile, ...
            'Resolution', 300, ...
            'BackgroundColor', 'current', ...
            'ContentType', 'image', ...   % avoid vector-crop quirks
            'Padding', 100);               % pixels of extra breathing room
        catch
            saveas(f, outFile);
        end
    end
end

% ====================================================================== %
% Helpers
% ====================================================================== %

function draw_swc_colored(ax, coord, id, pid, d, lw)
% Edges colored by endpoint distance (interpolated); nodes colored too.
    axes(ax); cla(ax,'reset'); hold(ax,'on');

    [~, pidx] = ismember(pid, id);
    e = find(pidx > 0);
    if ~isempty(e)
        j = pidx(e);
        % Each column is one segment (child->parent)
        X = [coord(e,1)'; coord(j,1)'];
        Y = [coord(e,2)'; coord(j,2)'];
        Z = [coord(e,3)'; coord(j,3)'];
        C = [d(e)';       d(j)'];

        h = surface(ax, X, Y, Z, C, ...
            'FaceColor','none', 'EdgeColor','interp', 'LineWidth', 2.*lw);
        % Prevent cross-column connections (avoids stray lines)
        if isprop(h,'MeshStyle'), set(h,'MeshStyle','column'); end
    end

    % Nodes colored by distance as well
    scatter3(ax, coord(:,1), coord(:,2), coord(:,3), 7, d, 'filled');

    grid(ax,'on'); set(ax,'YDir','normal');
end

function apply_graygrid(ax)
% Gray panel + light grids + white ticks/text (screenshot style).
    set(ax,'Color',          [0.40 0.40 0.40], ...
           'XColor',         'w', ...
           'YColor',         'w', ...
           'ZColor',         'w', ...
           'GridColor',      [0.85 0.85 0.85], ...
           'MinorGridColor', [0.65 0.65 0.65], ...
           'GridAlpha',      0.9, ...
           'MinorGridAlpha', 0.7, ...
           'TickDir',        'out', ...
           'Box',            'off', ...
           'Layer',          'top');
    grid(ax,'on'); grid(ax,'minor');
end

function match_right_axes_width(ax1, ax2)
% Force identical inner (drawable) widths and aligned left edges.
    set([ax1, ax2], 'Units','normalized', 'PositionConstraint','innerposition');
    ip1 = get(ax1,'InnerPosition'); ip2 = get(ax2,'InnerPosition');
    left  = max(ip1(1), ip2(1));   % align to farther-right left edge
    width = min(ip1(3), ip2(3));   % use the narrower inner width
    ip1(1) = left; ip1(3) = width;
    ip2(1) = left; ip2(3) = width;
    set(ax1,'InnerPosition',ip1); set(ax2,'InnerPosition',ip2);
end

function set_equal_cube_axes(axList, coord, padFrac)
% Apply the SAME cubic limits to all axes so scaling matches.
% padFrac adds symmetric padding (e.g., 0.05 = 5%).
    if nargin < 3, padFrac = 0.05; end
    mins = min(coord, [], 1);
    maxs = max(coord, [], 1);
    ctr  = (mins + maxs) / 2;
    span = max(maxs - mins);  if span == 0, span = 1; end
    span = span * (1 + padFrac*2);
    r = span/2;

    xl = [ctr(1)-r, ctr(1)+r];
    yl = [ctr(2)-r, ctr(2)+r];
    zl = [ctr(3)-r, ctr(3)+r];

    for ax = axList
        xlim(ax, xl); ylim(ax, yl); zlim(ax, zl);
    end
end