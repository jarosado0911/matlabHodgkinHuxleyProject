function P = hh_params(varargin)
%HH_PARAMS  Hodgkin–Huxley parameters (MKS units)
%   P = hh_params() returns defaults.
%   P = hh_params('gna',60e1,'gl',3e1) overrides fields by name-value.

  P = struct( ...
    'R'  , 250e-2 , ...   % ohm·m
    'C'  , 1e-2    , ...  % F/m^2  (≈1 µF/cm^2)
    'gk' , 5e1     , ...  % S/m^2
    'ek' , -90e-3  , ...  % V
    'gna', 50e1    , ...  % S/m^2
    'ena', 50e-3   , ...  % V
    'gl' , 0.0e1   , ...  % S/m^2
    'el' , -70e-3  ...    % V
  );

  % Optional overrides
  if mod(nargin,2) ~= 0
    error('Use name/value pairs, e.g., hh_params(''gl'',3e1).');
  end
  for k = 1:2:numel(varargin)
    name = varargin{k};  val = varargin{k+1};
    if ~isfield(P,name), error('Unknown parameter "%s".',name); end
    P.(name) = val;
  end
end
