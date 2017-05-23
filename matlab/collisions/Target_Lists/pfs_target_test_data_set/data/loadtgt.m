function tgt=loadtgt(filename, fieldnum)
% tgt = loadtgt(filename, fieldnum)
% loads field [0-9] from filename (modified ETS output)

    etstgts = loadxyz(filename,'obs','f','x','y','ra','dec');
    fibers = loadxyz('fiberpos.el90.pfi.dat','k','x','y','r');

    % positions and fiber indices of targets
    tgt_xy = etstgts.c(etstgts.obs == fieldnum);
    fibindx = etstgts.f(etstgts.obs == fieldnum);

    tgt = fibers.c + 1.5*i; % initialize targets at cobra centers
    tgt(fibindx) = tgt_xy; % transfer assigned targets to output
    tgt = conj(tgt*exp(i*pi/2)); % make transformation to simFun coordinates
    
    plot(conj(tgt)*exp(-i*pi/2) - fibers.c,'.');
    hold on;
    title([filename ': target location in Cobra patrol area'])
    axis equal;
    patrol_xy = conj(tgt)*exp(-i*pi/2) - fibers.c;
    % print('hint: q = simFun(tgt,''full'',1,0)')