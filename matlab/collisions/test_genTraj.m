%% test code to set up inputs for generateTrajectory

% bench definition
centers = getCentersRails(14);
centers = bsxfun(@times, [1 exp(i*2*pi/3) exp(i*4*pi/3)], centers);
centers = reshape(centers,[],1);
bench = defineBenchGeometry(centers,1,1);
bench.alpha = 0;
ncobras = length(centers);
clear centers

dens = .1:.1:5; % # of targets per patrol area
dist = 0:.2:5; % distances

for jj = 1:length(dens)
    tgt = assign_targets(dens(jj),bench);
    drawnow;
    dst(jj,:) = full(tgt.srtDIST(:,1)).';
end
dst(dst==0) = nan;
for jj = 1:length(dens)
    ntgts = sum(~isnan(dst(jj,:)));
    hh(jj,:)  = hist(dst(jj,:),dist)/ntgts;
end
figure(1)
imagesc(dist,dens,hh)
hold on;
plot(nanmean(dst,2), dens, 'k','linewidth',5)
hold off;



%% generate trajectory test code
% target definition
% $$$ THT = rand(2394,1)*2*pi;
% $$$ dA = 1 ./ ( (bench.rMax./bench.rMin).^2 - 1 ); 
% $$$ rRange = sqrt(bench.rMax.^2 - bench.rMin.^2);
% $$$ RDS = bsxfun(@times, sqrt(bsxfun(@plus, rand(size(THT)), dA)), rRange);
% $$$ targets   = bsxfun(@plus, RDS.*exp(1i*THT), bench.center);
% $$$ clear dA rRange THT RDS
% $$$ 
% $$$ gt2 = generateTrajectory2(targets(:,1), bench,1);


