function [l,c]=lines_to_skip(filename)
%filename = 'morphos/578-2-1-R-540-3530.CNG.swc';

fid = fopen(filename,'r');
l = 0;
tline = fgetl(fid);

% count how many lines to skip that either have '#' first or
% are empty lines
fprintf('Determining how many lines to skip in .swc file...\n')
while ischar(tline)
  
  if ((isempty(tline)) || (tline(1)=='#'))
    l = l+1;
  end  
  
  tline = fgetl(fid);
end
fclose(fid);
fprintf(sprintf('In file %s\n', filename));
fprintf(sprintf('skip the first %i lines...\n',l));

fid = fopen(filename,'r');
for i=1:l+1
    tline=fgetl(fid);
end

if (tline(1) == ' ') 
    c=1;
else
    c=0;
end

fclose(fid);

end

