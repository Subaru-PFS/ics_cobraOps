function output = applyHorns(data,refpos)

tref=[];
for ii=data.fixedCobras
    fldID = sprintf('pId%d',ii);
    tref = [tref data.(fldID).CCDpos];
end

hh = horn87(tref,refpos);

for ii=data.activeCobras
    fldID = sprintf('pId%d',ii);
    data.(fldID).CMSpos = data.(fldID).CCDpos.* exp(i*hh.R) .* hh.S + hh.T;
end

output = data;
end