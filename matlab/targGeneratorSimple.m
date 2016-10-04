function targGenerator(ntargets)
%TARGET GENERATOR
  % total targets added to file(s) will be ntargets.
if ~exist('ntargets','var'), ntargets = 1; end
 
   figure(1);
   hold on;
    bench = defineBenchGeometry([], 1,1);
    
            % angular coordinate of target
            numPos = length(bench.center); 
            THT = rand(numPos,ntargets)*2*pi;

            % radial coordinate of target
            dA = 1 ./ ( (bench.rMax./bench.rMin).^2 - 1 ); %fraction of the keepout area.
            rRange = sqrt(bench.rMax.^2 - bench.rMin.^2);
            RDS = bsxfun(@times, sqrt(bsxfun(@plus, rand(size(THT)), dA)), rRange);

            %% Rule 2 No targets closer than 2mm to any line in its target location.
             targets   = bsxfun(@plus, RDS.*exp(1i*THT), bench.center);
    keyboard;
    for jj = 1:length(bench.pids)
        fname = sprintf('TargetList_mId_%d_pId_%d.txt',bench.mids(jj),bench.pids(jj));
        tfile = fopen(fname,'a');
      

           %    targetoutput = [targets, zeros(length(targets)), zeros(length(targets))]; 
             plot(targets(jj,:),'b.');
            plot(bench.center(jj),'rx');
            bench.pids(jj)
            for qq = 1:length(targets(jj,:))
                fprintf(tfile,'%.2f,%.2f, 0, 0\n',real(targets(jj,qq)),imag(targets(jj,qq)));
            end
           %        fprintf(1,'wrote %s\n',fname)
        fclose(tfile);
    end
end

