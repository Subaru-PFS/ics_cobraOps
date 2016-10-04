function output=generateTargets(density,bench)
% usage output=generate_targets(density,bench)
% density is the # targets per patrol area
% specify bench if you have one, otherwise a full bench will be generated
%
% generate a uniform density set of targets over the field of view
% of the bench
%
% the following will generate an as-is bench from the CfgXml file:
% bench = defineBenchGeometry([],1,1)


toggle.makefigs = true;


% bench definition
if ~exist('bench','var')
    centers = getCentersRails(14);
    centers = bsxfun(@times, [1 exp(i*2*pi/3) exp(i*4*pi/3)], centers);
    centers = reshape(centers,[],1);
    bench = defineBenchGeometry(centers,1,1);
    bench.alpha = 0;
    clear centers
end
ncobras = length(bench.center);

%% uniformly generate targets over field of view.
ntgt = ceil(density * (bench.field.R/median(bench.rMax))^2);
% $$$ fprintf(1,'Generating %d targets for %d positioners\n',ntgt,ncobras);

THT = rand(1,ntgt) * 2 * pi;
RR  = sqrt(rand(1,ntgt)) * bench.field.R;
tgt = RR .* exp(i*THT); % generated targets, centered on cobra CM
tgt = tgt + bench.field.cm;

output = tgt;