function targGenerator(ndiff,nsame, usewait)
%TARGET GENERATOR
% ndiff = number of different targets
% nsame = number of the same target
% usewait = use wait moves for collision avoidance.
% total targets added to file(s) will be ndiff X nsame.
if ~exist('ndiff','var'), ndiff = 1; end
if ~exist('nsame','var'), nsame = 1; end
if ~exist('usewait','var'), usewait = 0; end

for dd = 1:ndiff
    simout = simFun(1,'',1,1,'alpha',0);
    bench = simout.bench;
    nSteps = simout.nSteps;
    targetoutput = [simple(simout.targets) round(nSteps.dtht), round(nSteps.dphi)];
    for jj = 1:length(bench.pids)
        fname = sprintf('TargetList_mId_%d_pId_%d.txt',bench.mids(jj),bench.pids(jj));
        tfile = fopen(fname,'a');
        for ss=1:nsame
            if(usewait)
                fprintf(tfile,'%.2f,%.2f,%d,%d\n',targetoutput(jj,:));
            else
                fprintf(tfile,'%.2f,%.2f, 0, 0\n',targetoutput(jj,1:2));
            end
        end
        %        fprintf(1,'wrote %s\n',fname)
        fclose(tfile);
    end
end

