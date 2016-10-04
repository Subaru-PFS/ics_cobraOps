function data = readsdf(filename)
% READSDF: read surface profile file from microXAM
% USAGE: data = readsdf(filename)
% 
% appends '.SDF' if not present

%append '.SDF' suffix if necessary
  strrep(filename,'.sdf','.SDF');
  if isempty(strfind(filename,'.SDF'))
    filename = [filename '.SDF'];
  end
  fid = fopen(filename);

  % READ HEADER
  skiprow = 1;
  fline = fgetl(fid);
  while (fline ~= '*')
    skiprow = skiprow + 1;
    %count words
    [start_word,end_word] = regexp(fline,'[^\s]+');
    if (length(start_word) > 2)
      value = fline(start_word(3):end_word(3));
      if ~isempty(str2num(value))
        %delete leading whitespace
        fline = regexprep(fline,'^\s+','');
        eval(['data.header.' fline ';']); 
      end
    end
    fline = fgetl(fid);
  end
  
  ncol = data.header.NumPoints;
  nrow = data.header.NumProfiles;
  ndata = nrow * ncol;
  % find num data points/line
  while isempty(fgetl(fid))
    skiprow = skiprow + 1;
  end
  fline = fgetl(fid);
  [start,finish] = regexp(fline,'[^\s]+');
  nperline = length(start);
  fclose(fid);

  %read in data
  data.header.skiprow = skiprow;
  data.header.nperline = nperline;
  R2 = ndata/nperline + skiprow - 1;
  C2 = nperline - 1;
  m = dlmread(filename,' ',[skiprow 0 R2 C2]);
  m = reshape(m', ncol, nrow)'*data.header.Zscale;
  % reshape array and make x,y vectors
  data.m = m(2:nrow-1,2:ncol-1);
  data.x = (0:ncol-3)*data.header.Xscale;
  data.y = (0:nrow-3)'*data.header.Yscale;
  return
  
  fclose(fid);
  