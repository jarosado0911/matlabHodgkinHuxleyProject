function makematlabmovtraces(dataFolder, geometryFile, outname, idx)
% Left: neuron scatter colored by per-node voltage (vm_t*.dat) + COLORED EDGES
% Right: four stacked time-series (rec_u, rec_h, rec_n, rec_m) for row idx
% Loads traces from trace_data.mat in dataFolder.

    if nargin < 4 || isempty(idx), idx = 1; end

    % ---------------- Geometry ----------------
    % need id & pid to build edges
    [~, id, pid, coords, r, ~]  = readswc(geometryFile);
    markerSize = (r./max(r)) * 50;

    % child->parent connectivity (indices into coords)
    [~, pidx] = ismember(pid, id);   % parent index for each node (0 if none)
    E = find(pidx > 0);              % children that have a parent
    J = pidx(E);                     % their parent indices

    % columns = segments (child->parent)
    X = [coords(E,1)'; coords(J,1)'];
    Y = [coords(E,2)'; coords(J,2)'];
    Z = [coords(E,3)'; coords(J,3)'];

    % ---------------- Load traces ----------------
    trFile = fullfile(dataFolder,'trace_data.mat');
    if ~isfile(trFile)
        error('Could not find %s. Save traces to trace_data.mat first.', trFile);
    end
    tr = load(trFile);  % expects t, rec_u, rec_h, rec_m, rec_n, record_index

    t     = tr.t(:);
    rec_u = tr.rec_u;    % B x T
    rec_h = tr.rec_h;    % B x T
    rec_n = tr.rec_n;    % B x T
    rec_m = tr.rec_m;    % B x T
    
    [~,T] = size(rec_u);
    if numel(t) ~= T
        error('Length of t (%d) does not match time dimension of rec_* (%d).', numel(t), T);
    end

    % ---------------- Figure & layout ----------------
    fig = figure('Units','normalized','Position',[0 0 1 1], ...
                 'Color',[0 0 0], 'InvertHardcopy','off');

    % layout knobs
    L=0.06; R=0.08; Bm=0.08; Tm=0.08; gapLR=0.03; gapY=0.04;
    leftW = 0.54; rightW = 1 - (L+R+leftW+gapLR);

    axLeft = axes('Parent',fig,'Units','normalized', ...
        'Position',[L, Bm, leftW, 1-(Bm+Tm)], 'PositionConstraint','innerposition');

    % Right column: 4 stacked axes
    hAvail = 1-(Bm+Tm) - 3*gapY;
    hEach  = hAvail/4;
    xRight = L + leftW + gapLR;

    axU = axes('Parent',fig,'Units','normalized','Position',[xRight, Bm + 3*(hEach+gapY), rightW, hEach], 'PositionConstraint','innerposition');
    axH = axes('Parent',fig,'Units','normalized','Position',[xRight, Bm + 2*(hEach+gapY), rightW, hEach], 'PositionConstraint','innerposition');
    axN = axes('Parent',fig,'Units','normalized','Position',[xRight, Bm + 1*(hEach+gapY), rightW, hEach], 'PositionConstraint','innerposition');
    axM = axes('Parent',fig,'Units','normalized','Position',[xRight, Bm + 0*(hEach+gapY), rightW, hEach], 'PositionConstraint','innerposition');

    % ---------------- Left plot (init with edges + nodes) ----------------
    cmin_mV = -5;  cmax_mV = 50;
    u0 = readmatrix(fullfile(dataFolder, 'data', sprintf('vm_t%d.dat', 0)));
    u0 = u0(:);

    axes(axLeft); cla(axLeft); hold(axLeft,'on');

    % edges colored by endpoint values (interpolated along segment)
    C0 = [u0(E)'; u0(J)'];  % 2 x numSegments
    hEdges = surface(axLeft, X, Y, Z, C0, ...
        'FaceColor','none', 'EdgeColor','interp', 'LineWidth', 2.0);
    if isprop(hEdges,'MeshStyle'), set(hEdges,'MeshStyle','column'); end

    % nodes
    hDots = scatter3(axLeft, coords(:,1), coords(:,2), coords(:,3), ...
                     markerSize, 'filled', 'CData', u0);

    focusNode = idx;                  % coordinate index to highlight
    focusSize = max(120, 4*markerSize(focusNode));% make it stand out
    hFocus = scatter3(axLeft, ...
        coords(focusNode,1), coords(focusNode,2), coords(focusNode,3), ...
        focusSize, 'o', 'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor','w', 'LineWidth',1.5);
    uistack(hFocus,'top');   

    view(axLeft, 2);  % XY
    colormap(fig, 'jet');
    caxis(axLeft, [cmin_mV cmax_mV]*1e-3);
    title(axLeft, sprintf('t = %0.2f ms', t(1)*1e3), 'Color','w','FontWeight','bold');
    xlabel(axLeft, '\mu m', 'Color','w'); ylabel(axLeft, '\mu m', 'Color','w');
    style_graypanel(axLeft);

    % ---- horizontal colorbar INSIDE top-right of axLeft ----
    drawnow;
    cb = colorbar(axLeft);
    cb.Orientation   = 'horizontal';
    cb.Color         = 'w';
    cb.Label.String  = 'Voltage [mV]';
    cb.Label.Color   = 'w';
    cb.TickDirection = 'out';
    cb.Units         = 'normalized';
    cb.Location      = 'manual';

    axIP   = get(axLeft,'InnerPosition');
    pad    = 0.02;
    barLen = 0.35 * axIP(3);
    barThk = 0.015 * axIP(4);
    leftCB = axIP(1) + axIP(3) - barLen - pad;   % top-right inset
    botCB  = axIP(2) + axIP(4) - barThk - 3*pad;
    cb.Position = [leftCB, botCB, barLen, barThk];

    % ---------------- Right: time-series plots (init once) ----------------
    yU = 1e3 *  rec_u(1,:).';
    yH =        rec_h(1,:).';
    yN =        rec_n(1,:).';
    yM =        rec_m(1,:).';

    hold(axU,'on'); pU = plot(axU, t, yU, 'LineWidth', 1.5); pU.Color=[0.9 0.95 1]; style_graypanel(axU);
    ylabel(axU, sprintf('u @ node %d [mV]', idx), 'Color','w'); title(axU,'rec\_u','Color','w');
    xlim(axU, [t(1) t(end)]);

    hold(axH,'on'); pH = plot(axH, t, yH, 'LineWidth', 1.5); pH.Color=[1 0.9 0.6]; style_graypanel(axH);
    ylabel(axH, 'h [-]', 'Color','w'); title(axH,'rec\_h','Color','w'); ylim(axH,[0 1]); xlim(axH,[t(1) t(end)]);

    hold(axN,'on'); pN = plot(axN, t, yN, 'LineWidth', 1.5); pN.Color=[0.7 1 0.7]; style_graypanel(axN);
    ylabel(axN, 'n [-]', 'Color','w'); title(axN,'rec\_n','Color','w'); ylim(axN,[0 1]); xlim(axN,[t(1) t(end)]);

    hold(axM,'on'); pM = plot(axM, t, yM, 'LineWidth', 1.5); pM.Color=[0.8 0.8 1.0]; style_graypanel(axM);
    ylabel(axM, 'm [-]', 'Color','w'); xlabel(axM,'Time [s]','Color','w'); title(axM,'rec\_m','Color','w'); ylim(axM,[0 1]); xlim(axM,[t(1) t(end)]);

    % moving time cursor (xline) on each subplot
    xu = xline(axU, t(1), 'w-'); xh = xline(axH, t(1), 'w-'); xn = xline(axN, t(1), 'w-'); xm = xline(axM, t(1), 'w-');
    hotpink = [1, 0.4118, 0.7059]; set([xu xh xn xm], 'Color', hotpink, 'LineWidth', 2.5, 'LineStyle','-');

    % ---------------- Video writer ----------------
    v = VideoWriter(sprintf('%s.mp4', outname), 'MPEG-4');
    open(v);

    % ---------------- Animate ----------------
    for k = 1:2:T
        fileIdx = k-1; % files are vm_t0, vm_t1, ...
        u_sol = readmatrix(fullfile(dataFolder, 'data', sprintf('vm_t%d.dat', fileIdx)));
        u_sol = u_sol(:);

        % update node colors and edge colors
        set(hDots,  'CData', u_sol);
        set(hEdges, 'CData', [u_sol(E)'; u_sol(J)']);

        % update title and time cursors
        title(axLeft, sprintf('t = %0.2f ms', t(k)*1e3), 'Color','w','FontWeight','bold');
        xu.Value = t(k); xh.Value = t(k); xn.Value = t(k); xm.Value = t(k);

        drawnow limitrate
        writeVideo(v, getframe(fig));
    end

    close(v);
end

% ---- Styling helper: gray panel with gridlines and white ticks ----
function style_graypanel(ax)
    set(ax, 'Color', [0.40 0.40 0.40], ...
            'XColor','w','YColor','w','ZColor','w', ...
            'GridColor',[0.5 0.5 0.5], 'GridAlpha',0.9, ...
            'TickDir','out', 'Box','off', 'Layer','top', ...
            'XMinorGrid','off','YMinorGrid','off','ZMinorGrid','off');
    grid(ax,'on');
end