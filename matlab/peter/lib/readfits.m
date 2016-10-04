function fits=readfits(filename)
% a better wrapper for fitsinfo and fitsread?

info = fitsinfo(filename);

for jj=1:length(info.PrimaryData)
  kws = info.PrimaryData(jj).Keywords;
  for kk=1:length(kws)
    fldname = strrep(kws{kk,1},'-','_');
    info.PrimaryData(jj).kw.(fldname) = kws{kk,2};
  end
end

pdata = fitsread(filename);

if isfield(info,'BinaryTable')
  for jj=1:length(info.BinaryTable)
    kws = info.BinaryTable(jj).Keywords;
    for kk=1:length(kws)
      fldname = strrep(kws{kk,1},'-','_');
      info.BinaryTable(jj).kw.(fldname) = kws{kk,2};
    end

    fdata{jj} = fitsread(filename,'BinaryTable',jj);
    if ~isempty(fdata{jj})
      btname = ['b' num2str(jj)];
      for kk=1:length(kws)
        if regexp(kws{kk,1},'TTYPE')
          fname = lower(kws{kk,2});
          matches = regexp(kws{kk,1},'\d+','match');
          col = str2num(matches{1});
          if iscell(fdata{jj}{col})
            fdata{jj}{col} = strtrim(fdata{jj}{col});
          end
          fits.(btname).(fname) = fdata{jj}{col};
        end
      end
    end
  end
end
  
fits.info = info;
fits.data = pdata;