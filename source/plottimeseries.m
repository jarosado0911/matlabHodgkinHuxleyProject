function plottimeseries(t, v, i)
%PLOTTIMESERIES Plot the i-th trace of v versus time t in a new window.
%   plottimeseries(t, v, i)
%     - t : 1xT or Tx1 time vector (seconds)
%     - v : BxT or TxB matrix of series (e.g., u_rec, h_rec, m_rec, n_rec)
%     - i : row/series index to plot (1-based)
%
%   The function auto-detects whether time runs across columns or rows.

    % ---- validate t ----
    if ~isnumeric(t) || ~isvector(t)
        error('t must be a numeric vector.');
    end
    t = t(:);  % column vector
    T = numel(t);

    % ---- validate v orientation & index ----
    [r, c] = size(v);
    if c == T
        % time along columns: series are rows
        B = r;
        if ~(isscalar(i) && i>=1 && i<=B && i==floor(i))
            error('Index i must be an integer in [1, %d].', B);
        end
        y = v(i, :).';
        seriesCount = B;
        orientation = 'rows are series';
    elseif r == T
        % time along rows: series are columns
        B = c;
        if ~(isscalar(i) && i>=1 && i<=B && i==floor(i))
            error('Index i must be an integer in [1, %d].', B);
        end
        y = v(:, i);
        seriesCount = B;
        orientation = 'columns are series';
    else
        error('Size(v) = [%d x %d] is incompatible with numel(t) = %d.', r, c, T);
    end

    % ---- figure & axes with gray panel and gridlines ----
    f  = figure('Color',[0 0 0], 'InvertHardcopy','off'); %#ok<NASGU>
    ax = axes('Parent', gcf, 'Units','normalized', 'Position',[0.12 0.14 0.82 0.80]);
    set(ax,'Color',[0.40 0.40 0.40], ...
           'XColor','w','YColor','w', ...
           'GridColor',[0.85 0.85 0.85], 'MinorGridColor',[0.65 0.65 0.65], ...
           'GridAlpha',0.9, 'MinorGridAlpha',0.7, ...
           'TickDir','out', 'Box','off', 'Layer','top');
    grid(ax,'on'); grid(ax,'minor'); hold(ax,'on');

    % ---- plot ----
    p = plot(ax, t, y, 'LineWidth', 1.8);
    set(p, 'Color', [0.85 0.95 1.00]);  % soft cyan line
    axis(ax,'tight');

    % ---- labels & title ----
    vn = inputname(2);  % try to show the variable name (u_rec, h_rec, etc.)
    if isempty(vn), vn = 'series'; end
    xlabel(ax, 'Time [s]', 'Color','w');
    ylabel(ax, infer_ylabel(vn), 'Color','w');
    title(ax, sprintf('%s at index i = %d (%s, total %d)', vn, i, orientation, seriesCount), ...
          'Color','w','FontWeight','bold');

end

% --- helper: pick a nice y-label based on the variable name ---
function yl = infer_ylabel(vname)
    vn = lower(vname);
    if contains(vn,'u') && ~contains(vn,'hu') && ~contains(vn,'mu') && ~contains(vn,'nu')
        yl = 'Voltage [V]';
    elseif contains(vn,'h') || contains(vn,'m') || contains(vn,'n')
        yl = 'Gating variable [-]';
    else
        yl = 'Value';
    end
end
