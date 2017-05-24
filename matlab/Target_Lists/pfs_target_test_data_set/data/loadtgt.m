function tgt=loadtgt(filename, fieldnum, default_offset)
% tgt = loadtgt(filename, fieldnum, [default_offset=1.5i])
% loads field [0-9] from filename (modified ETS output)
% use default_offset=nan to generate figures.
% output is prepared to be input into simFun

    etstgts = load_ets_txt(filename);
    fibers = loadxyz('fiberpos.el90.pfi.dat','k','x','y','r');

    % positions and fiber indices of targets
    tgt_xy = etstgts.c(etstgts.obs == fieldnum);
    fibindx = etstgts.fib(etstgts.obs == fieldnum);

    % ra dec positions with ra corrected for declination.
    tgt_ra  = etstgts.ra(etstgts.obs == fieldnum);
    tgt_dec = etstgts.dec(etstgts.obs == fieldnum);
    tgt_radec = (tgt_ra - nanmean(tgt_ra)) .* cos(tgt_dec * pi/180) ...
        + nanmean(tgt_ra) + i*tgt_dec;

     % use 1.5*i for target lists use nan for distribution analysis
    if ~exist('default_offset','var'), default_offset = 1.5*i; end;

    tgt = fibers.c + default_offset; % initialize targets at cobra centers
    tgt(fibindx) = tgt_xy; % transfer assigned targets to output
    tgt = conj(tgt*exp(i*pi/2)); % make transformation to simFun coordinates
    
    patrol_xy = conj(tgt)*exp(-i*pi/2) - fibers.c;

    if isnan(default_offset)

        figure(1)
        plot(patrol_xy,'.');
        title(sprintf('%s[%d]: target location in Cobra patrol area',filename,fieldnum),'interpreter','none');
        xlabel('[mm]');
        ylabel('[mm]');
        axis equal;
        % print('hint: q = simFun(tgt,''full'',1,0)')

        figure(2)
        plot(tgt_radec,'.');
        hold on;
        title(sprintf('%s[%d]: target location on "sky"',filename,fieldnum),'interpreter','none');
        xlabel('RA "x cos(Dec)" [deg]');
        ylabel('Dec [deg]');
        axis equal;

        figure(3)
        histogram(abs(patrol_xy));
        title(sprintf('%s[%d]: target radial distribution in PA',filename,fieldnum),'interpreter','none');
        xlabel('r [mm]');
        ylabel('# fibers');
        hold on;
    end
    
end
    
function output=load_ets_txt(filename)
% loads ETS output file

% $$$   data format is:
% $$$   1  ID
% $$$   2  Fiber
% $$$   3  X
% $$$   4  Y
% $$$   5  R.A.          [deg.]
% $$$   6  Dec.          [deg.]
% $$$   7  observation #


   data_format = '%s %d %f %f %f %f %d';
   fid = fopen(filename);
   contents = textscan(fid, data_format,'CommentStyle','#');
   fclose(fid);
   
   output.obs = contents{7};
   output.fib = contents{2};
   output.x   = contents{3};
   output.y   = contents{4};
   output.ra  = contents{5};
   output.dec = contents{6};
   output.c   = output.x + i*output.y;
end