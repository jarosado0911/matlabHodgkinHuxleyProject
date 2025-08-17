function print_sim_params(S)
%PRINT_SIM_PARAMS  Pretty-print simulation parameter struct.

  Units = struct( ...
    'dt'     , 's', ...
    'endTime', 's', ...
    'delay'  , 's', ...
    'vStart' , 'V', ...
    'vClamp' , 'V', ...
    'nT'     , 'steps', ...
    'ni'     , '-', ...
    'mi'     , '-', ...
    'hi'     , '-' );

  fn = fieldnames(S);
  w  = max(cellfun(@numel, fn));
  fprintf('\nSimulation parameters:\n');
  for i = 1:numel(fn)
    k   = fn{i};
    val = S.(k);
    u   = Units.(k);
    if isscalar(val)
      if isnumeric(val)
        fprintf('  %-*s : %-12.6g %s\n', w, k, val, u);
      else
        fprintf('  %-*s : %s %s\n', w, k, mat2str(val), u);
      end
    else
      fprintf('  %-*s : %s %s\n', w, k, mat2str(val), u);
    end
  end
end
