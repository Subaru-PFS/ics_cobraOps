function output=readrawcats(cats)

% analyze position stability

  cat = ls2cell(cats);

  for jj=1:length(cat)
    data(jj) = loadxyz(cat{jj},'x','y','xw','yw','xmin','ymin','xmax','ymax','flux','thresh','max','area','fwhm');
  end

  output = data;
  end
