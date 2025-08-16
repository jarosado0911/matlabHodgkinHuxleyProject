function S = sim_params(dt, varargin)
%SIM_PARAMS  Simulation parameters for neuron runs (MKS units)
%   S = SIM_PARAMS(dt) returns defaults using the provided time step dt.
%   S = SIM_PARAMS(dt,'endTime',80e-3, 'vStart',-65e-3, ...) overrides fields.
%
% Fields:
%   dt, endTime, delay, vStart, vClamp, nT, ni, mi, hi

  if nargin < 1
    error('sim_params:MissingArg','Provide dt (seconds).');
  end

  % Defaults (matching your current code)
  S = struct( ...
    'dt'     , dt, ...        % s
    'endTime', 50e-3, ...     % s
    'delay'  , 15e-3, ...     % s
    'vStart' , 0e-3, ...      % V
    'vClamp' , 50e-3, ...     % V
    'ni'     , 0.0376969, ... % gating initial (dimensionless)
    'mi'     , 0.0147567, ... % gating initial (dimensionless)
    'hi'     , 0.995941 ...   % gating initial (dimensionless)
  );

  % Overrides (name/value)
  if mod(numel(varargin),2) ~= 0
    error('sim_params:NameValue','Use name/value pairs for overrides.');
  end
  for k = 1:2:numel(varargin)
    name = varargin{k};
    val  = varargin{k+1};
    if ~isfield(S,name)
      error('sim_params:UnknownField','Unknown field "%s".',name);
    end
    S.(name) = val;
  end

  % Derived
  S.nT = floor(S.endTime / S.dt);  % number of time steps
end
