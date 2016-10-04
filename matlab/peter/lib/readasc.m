function data = readasc(filename)
% READASC: read single line surface profile file from microXAM
% USAGE: data = readasc(filename)
% data is a structure with fields x and y

  %add '.asc' suffix if necessary
  if isempty(strfind(filename,'.asc'))
    filename = [filename '.asc'];
  end
  fid = fopen(filename);

  % MEASURE HEADER
  skiprows = 0;
  fline = fgetl(fid);
  while ~isempty(regexp(fline,'^#','once'))
    skiprows = skiprows + 1;
    fline = fgetl(fid);
  end
  fclose(fid);
  M = dlmread(filename,'\t',skiprows,0);
  nrows = length(M);
  data.x.data = M(:,1) * nrows/(nrows+1); % adjust spacing for their error
  data.x.unit = '\mum';
  data.x.label = 'profile coordinate';
  data.y.data = M(:,2)*1e-3; % microns
  data.y.unit = '\mum';
  data.y.label = 'surface height';
  return
