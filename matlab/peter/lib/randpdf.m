function output=randpdf(pdf,nevents)
% generates a random distribution based on a probability
% distribution specified by PDF

pdfsize = size(pdf);

vlength = numel(pdf);
pdfvec = reshape(pdf, 1, vlength);

cdf = cumsum(pdfvec) / sum(pdfvec);

if matlabpool('size') > 0
  parfor jj=1:nevents
    rand_dist(jj) = find(rand < cdf, 1) - 1;
  end
else
  for jj=1:nevents
    rand_dist(jj) = find(rand < cdf, 1) - 1;
  end
end
rdf = hist(rand_dist,1:vlength);
output = reshape(rdf, pdfsize);
