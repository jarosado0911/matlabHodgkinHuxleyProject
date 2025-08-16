function print_hh_params(P)
%PRINT_HH_PARAMS  Pretty-print HH parameter struct to console.

  % Optional: units for known fields
  Units = struct( ...
      'R'  , 'ohm·m', ...
      'C'  , 'F/m^2', ...
      'gk' , 'S/m^2', ...
      'ek' , 'V', ...
      'gna', 'S/m^2', ...
      'ena', 'V', ...
      'gl' , 'S/m^2', ...
      'el' , 'V' );

  f = fieldnames(P);
  w = max(cellfun(@numel, f));       % column width for names

  fprintf('\nHodgkin–Huxley parameters:\n');
  for k = 1:numel(f)
      name = f{k};
      val  = P.(name);
      unit = '';
      if isfield(Units, name)
          unit = Units.(name);
      end
      if isscalar(val)
          fprintf('  %-*s : %-12.6g %s\n', w, name, val, unit);
      else
          fprintf('  %-*s : %s %s\n', w, name, mat2str(val), unit);
      end
  end
end
