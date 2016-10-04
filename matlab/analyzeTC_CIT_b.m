function output=analyzeTC_CIT_b(data,cal,maxsteps,beta)
% analyze target convergence with manual entry of beta
% output=analyzeTC_CIT_b(data,cal,maxsteps,beta)
% data: *_str data structure from Chaz
% cal : pId# structure from ../../../SVN/MATLAB/EMconfig_MM-DD-YY.mat
% maxsteps: ignore convergences taking more steps than maxsteps
%
% output: J1, J2 error/request arrays,
% if REDUCE works, J1r and J2r hold the alpha and beta values for
% the motors.
%
% if the motor is underdamped, there is commented out code to remove
% the damping with a line fit to error vs. request.

  if ~exist('maxsteps','var'), maxsteps = 16; end
  
  %% motor calibration constants
  DtoR = pi/180;
  PDE = 1.0; %um
  pix2um = 90; % ugh! hard coded...
  theta0 = str2num(cal.KINEMATICS.Global_base_ori_z.Text) * DtoR;
  cobracenter = str2num(cal.KINEMATICS.Global_base_pos_x.Text) + str2num(cal.KINEMATICS.Global_base_pos_y.Text)*i;
  arm1 = str2num(cal.KINEMATICS.Link1_Link_Length.Text) * pix2um;
  arm2 = str2num(cal.KINEMATICS.Link2_Link_Length.Text) * pix2um;
  RR = @(phi) sqrt(cos(phi) * 2 * arm1 * arm2 + arm1^2 + arm2^2);

  ntargets = length(data);
  kk = 1;
  %% init outputs to column vectors
  J1.pos = [0;0];
  J1.req = [0;0];  % all moves request
  J1.err = [0;0];  % all moves error
% $$$   J1.req1 = [0;0]; % first move request
% $$$   J1.err1 = [0;0]; % first move error
  J1.req_pde = [0;0]; % angular error
  J1.err_pde = [0;0]; % angular error
  J1.vec = [0 0 ; 0 0 ];
  J2 = J1;
  BETA = 1; % 1 for johannes, .5 for peter
  
  for jj=1:ntargets
    if 1 %data(jj).cvrgIter <= maxsteps
      stat0 = [find(data(jj).status == 0) ;...
               find(data(jj).status == 1, 1)];
      kkstrt = stat0(1:end-1);
      kkfnsh = stat0(2:end);
      nsteps = length(kkstrt);
      kkrange = kk:(kk+nsteps-1);
      Rpos   = RR(data(jj).J2(stat0)); % fiber radial position in
                                       % patrol region
      %% throw the baby out with the bath water
      if (max(abs([data(jj).J1err(kkfnsh) ; data(jj).J2err(kkfnsh)])) < 1.5) 
        
        %% THETA MOTOR
        J1ang_i = data(jj).J1(kkstrt);
        J1ang_f = data(jj).J1(kkfnsh);
        J1req = data(jj).J1err(kkstrt);
        J1err = data(jj).J1err(kkfnsh);
        %% req vs. pos vectors
        strt = J1ang_i + i * J1err./abs(J1req).^BETA;
        fnsh = J1ang_f + i * J1err./abs(J1req).^BETA;
        J1.vec(kkrange,:) = [NaN NaN; strt(2:end) fnsh(2:end)];
        %% reqest and error angles
        J1.req(kkrange) = J1req;
        J1.err(kkrange) = J1err;
% $$$         %% reqest and error angles -- first moves
% $$$         J1.req1(jj) = data(jj).J1err(1);
% $$$         J1.err1(jj) = data(jj).J1err(2);
        %% angular uncertainties for fit weighting
        J1.req_pde(kkrange) = PDE./Rpos(kkstrt);
        J1.err_pde(kkrange) = PDE./Rpos(kkfnsh);
        
        %% PHI MOTOR
        J2ang_i = data(jj).J2(kkstrt);
        J2ang_f = data(jj).J2(kkfnsh);
        J2req = data(jj).J2err(kkstrt);
        J2err = data(jj).J2err(kkfnsh);
        %% req vs. pos vectors
        strt = J2ang_i + i * J2err./abs(J2req).^BETA;
        fnsh = J2ang_f + i * J2err./abs(J2req).^BETA;
        J2.vec(kkrange,:) = [NaN NaN; strt(2:end) fnsh(2:end)];
        %% reqest and error angles
        J2.req(kkrange) = J2req;
        J2.err(kkrange) = J2err;
% $$$         %% reqest and error angles -- first moves
% $$$         J2.req1(jj) = data(jj).J2err(1);
% $$$         J2.err1(jj) = data(jj).J2err(2);
        %% angular uncertainties for fit weighting
        J2.req_pde(kkrange) = PDE./Rpos(kkstrt);
        J2.err_pde(kkrange) = PDE./Rpos(kkfnsh);
        
        kk = kk + nsteps;
      end 
    end
  end
  
  okp1 = find(~isnan(J1.vec(:,1)) & J1.req > 0);
  okn1 = find(~isnan(J1.vec(:,1)) & J1.req < 0);
  okp2 = find(~isnan(J2.vec(:,1)) & J2.req > 0);
  okn2 = find(~isnan(J2.vec(:,1)) & J2.req < 0);
  
  eva.p1 = J1.vec(okp1,:).';
  eva.n1 = J1.vec(okn1,:).';
  eva.p2 = J2.vec(okp2,:).';
  eva.n2 = J2.vec(okn2,:).';
  
  J1 = modify_damping(J1);
  J2 = modify_damping(J2);

  try
    J1r = reduce(J1,[],beta);
    J2r = reduce(J2,[],beta);
    clear output;
    output = packstruct(J1r, J2r, J1, J2, eva);

    fh1 = figure('Name',sprintf('PID%s Phi',cal.DATA_HEADER.Positioner_Id.Text));
    subplot(2,3, [1 2]);
    make_loglogscatter(J2r,'phi','');
    make_scaledhist(J2r,'phi','');
    saveas(gcf,sprintf('PID%s_PhiModel',cal.DATA_HEADER.Positioner_Id.Text),'fig');
    
    fh2 = figure('Name',sprintf('PID%s Theta',cal.DATA_HEADER.Positioner_Id.Text));
    subplot(2,3, [1 2]);
    make_loglogscatter(J1r,'theta');
    make_scaledhist(J1r,'theta');
    saveas(gcf,sprintf('PID%s_ThetaModel',cal.DATA_HEADER.Positioner_Id.Text),'fig');
    
  catch err
    disp(err.message)
    output = packstruct(J1, J2, eva);
  end
end

%% correct damping factor
function Jn = modify_damping(Jn)
  Jn.err0 = Jn.err;
  pp = polyfit(Jn.req,Jn.err0,1); % line fit to get starting point
  % fit slope with y-int == 0;
  fitslope = fit(Jn.req, Jn.err0, @(a,x) a*x, 'Startpoint', pp(1));
  Jn.err = Jn.err0 - feval(fitslope,Jn.req);
  Jn.slope = fitslope;
  
  %% pos and neg slopes
  pos = find(Jn.req > 0);
  neg = find(Jn.req < 0);
  posslope = fit(Jn.req(pos), Jn.err0(pos), @(a,x) a*x, 'Startpoint', pp(1));
  if length(neg)>1
    negslope = fit(Jn.req(neg), Jn.err0(neg), @(a,x) a*x, 'Startpoint', pp(1));
  else
    negslope.a = NaN;
  end
  
  Jn.posslope = posslope;
  Jn.negslope = negslope;
  
end


%% REDUCE(data,cut)
function dataout = reduce(data,cut,beta)
% data is in the form of a signed req/err set
  if ~exist('cut','var') || isempty(cut)
    cut = logical(ones(size(data.req)));
  end
  
  %% for the envelope (beta), it is proper to fit the data in log space
  %  because we are interested in the ratios of data to model vs. request angle.
  p1_weights = sqrt( log10(data.err_pde(cut)./abs(data.err(cut)) + 1).^2 + ...
                     log10(data.req_pde(cut)./abs(data.req(cut)) + 1).^2).^(-1);

  fit_p1 = fit(log10(abs(data.req(cut))), log10(abs(data.err(cut))), 'poly1',...
               'weights',p1_weights);

  if ~exist('beta','var'), beta = fit_p1.p1; end;
  
  err_scaled = data.err(cut) ./ abs(data.req(cut)).^beta;

  histx = linspace(0, std(err_scaled) * 6, 100)'; %(0:.01:2)';
  h_err_scaled = histc( abs(err_scaled), histx );
  x_err_scaled = histx + 0.5*(histx(2)-histx(1));
  try
    roi = subarray(find(h_err_scaled > 1), 1:40);
  catch
    roi = find(h_err_scaled > 1);
  end

  %% for the scaled-error distribution function, it is better to fit the data
  %  in linear space while avoiding domains with sparse data 

  fit_e1 = fit(x_err_scaled(roi), h_err_scaled(roi), 'exp1', ...
               'weights', 1./h_err_scaled(roi));
  gauss1 = fittype(@(A,s,x) A/(s*sqrt(2*pi)).*exp(-0.5*(x./s).^2));
  fit_g1 = fit(x_err_scaled(roi), h_err_scaled(roi), gauss1,...
               'weights', 1./h_err_scaled(roi), 'Start', [1 1]);

  dataout.req             = data.req(cut);
  dataout.err             = data.err(cut);
  dataout.req_pde         = data.req_pde(cut);
  dataout.err_pde         = data.err_pde(cut);
  dataout.err_scaled      = err_scaled;
  dataout.hist.dens       = h_err_scaled;
  dataout.hist.err_scaled = x_err_scaled;
  dataout.alpha           = -1./fit_e1.b;
  dataout.beta            = beta;
  dataout.fitline         = fit_p1;
  dataout.fitexp1         = fit_e1;
  dataout.fgauss1         = fit_g1;

end

%% DATA=MAKEFAKERAWDATA(REQ, ALPHA, BETA, NOISE_FLOOR)
function data=MakeFakeRawData(req, alpha, beta, noise_floor)
% generate fake data from alpha, beta parameters
% this is the generalized case
  if ~exist('noise_floor','var'), noise_floor = 0; end;
  data.req = req;
  data.err = -alpha * log(1 - rand(size(req))) .* abs(req).^beta + ...
      randn(size(req)) * noise_floor;
  data.err = abs(data.err) .* sign(rand(size(req)) - .5);
end

% $$$ function fakedata=makefakedata(data, newbeta, req_xover)
% $$$ % make fake data with new beta based on existing data
% $$$ % preserve the crossover point
% $$$   if ~exist('req_xover','var')
% $$$     req_xover = (4.6*data.alpha)^(1/(1-data.beta));
% $$$   end
% $$$   req = data.req;
% $$$   newalpha = data.alpha * req_xover^data.beta / req_xover^newbeta;
% $$$   fakedata = MakeFakeRawData(req, newalpha, newbeta);
% $$$ end

%% MAKE_LOGLOGSCATTER(DATA,MOTORNAME,SUBTITLE)
%% esitmate erro bars for this one.
function make_loglogscatter(data,motorname,subtitle)

  if ~exist('motorname','var'), motorname = 'xi'; end;
  if ~exist('subtitle','var'), subtitle = ''; end;
  p(1) = errorbars(@loglog, ...
                   abs(data.req), data.req_pde, ...
                   abs(data.err), data.err_pde, 'b');
  hold on;

  xlimits = xlim;
  xx = logspace(log10(xlimits(1)), log10(xlimits(2)));
  p(2) = plot(xx, 10.^feval(data.fitline, log10(xx)),'r','LineWidth',2);
  
  p(3) = plot(xx,xx,'g','LineWidth',2);
% $$$   refline(0.1,0,Inf,'k--'); 
  hold off;
  xlabel(['|\' motorname '_{req}| [rad]']);%,'fontsize',14);
  ylabel(['|\Delta\' motorname '| [rad]']);%,'fontsize',14);
  lh = legend(p, 'Cobra moves', ...
              sprintf('\\Delta\\%s \\propto \\%s^{%.3f}',...
                      motorname,motorname,data.beta),...
              'break-even line','location','NW');
  set(lh,'fontsize',14);
  title(['Cobra \' motorname '  motion' subtitle],'fontsize',14);
  axis([1e-5 10 1e-5 10]);
end

%% MAKE_SCALEDHIST(DATA, MOTORNAME, SUBTITLE)
function make_scaledhist(data, motorname, subtitle)
  
  if ~exist('motorname','var'), motorname = 'xi'; end;
  if ~exist('subtitle','var'), subtitle = ''; end;

  set(0,'DefaultAxesFontSize','factory');%14);
  set(0,'DefaultTextFontSize','factory');%14);
  
  wLength = 100;
  ww = ones(wLength,1)/wLength;
  [xx, indx] = sort(abs(data.req));
  yy = data.err_scaled(indx);
  yy2 = yy.^2;
  xxc = conv(xx,ww,'valid');
  yyc = conv(yy,ww,'valid');
  yy2c = conv(yy2,ww,'valid');
  
  yy_std = sqrt( yy2c -yyc.^2 );
  yy_eim = yy_std / sqrt(wLength - 1);

  %% LOWER LEFT PLOT
  subplot(2,3,[4 5])
  h(1) = semilogx(abs(data.req), data.err_scaled, '.');  %smartzoom(5);
  aax = axis;
  hold on;
  errorbars(@semilogx, abs(data.req), data.req_pde,...
            data.err_scaled, data.err_pde ./ (abs(data.req).^data.beta));
  h(2) = plot(xxc,yyc,'r','LineWidth',2);
         plot(xxc,yyc+1*yy_std,'g','LineWidth',2);
  h(3) = plot(xxc,yyc-1*yy_std,'g','LineWidth',2);
 
% $$$          plot(xxc,yyc+3*yy_eim,'m','LineWidth',2);
% $$$   h(4) = plot(xxc,yyc-3*yy_eim,'m','LineWidth',2);
  hold off;
% $$$   legend(h, 'data','running mean','running std. dev. (1-\sigma)', ...
% $$$          'error in running mean (3-\sigma)','location','best');
   xlabel(['|\' motorname '_{req}| [rad]']); 
  ylabel(sprintf('\\%s_{err}/\\%s_{req}^{\\beta}',motorname, motorname));
  title('Scaled error vs. request angle');
  axis(aax);
  p1yrange = ylim;
  grid on;
  histx.ds = linspace(p1yrange(1),p1yrange(2)  ,1e2);
  hdata.ds = hist( data.err_scaled, histx.ds );

  %% LOWER RIGHT PLOT
  subplot(236)
  barh(histx.ds,hdata.ds/(sum(hdata.ds)*(histx.ds(2)-histx.ds(1))),'hist');
  ylim(p1yrange);
  xlabel('Prob. dens. (per unit scaled error)');
  ylabel(sprintf('\\%s_{err}/\\%s_{req}^{\\beta}',motorname,motorname)); 
  title('Scaled error distribution (linear)');
  grid on;
  
  %% UPPER RIGHT PLOT
  subplot(233)
  bar(data.hist.err_scaled, log(data.hist.dens),'BaseVal',-1);
  ylimit = ylim; 
  hold on;
  XX = data.hist.err_scaled;
  plot(XX, log(feval(data.fitexp1,XX)), 'r','LineWidth',2);
  plot(XX, log(feval(data.fgauss1,XX)), 'y','LineWidth',2);
  plot(XX, log(feval(data.fgauss1,XX)), 'k:','LineWidth',2);
  hold off;

  xlim([0 8*data.alpha]);
  ylim([-1 ylimit(2)]);

% $$$   legend('|error distribution|', ...
% $$$          sprintf('exponential fit, \\alpha = %.2f', -1/data.fitexp1.b), ...
% $$$          sprintf('gaussian fit, \\sigma = %.2f', data.fgauss1.s ));
  xlabel(sprintf('X = |\\%s_{err}/\\%s_{req}^{\\beta}|',motorname,motorname)); 
  ylabel('ln(#)');
  title('Scaled error distribution (ln vs abs)');
  t = text(.5,.25,sprintf('$\\propto \\exp(-\\frac{X}{%.2f})$', -1/data.fitexp1.b),...
           'Color','r','fontsize',20,'horizontal','left','units','normal');
  set(t,'Interpreter','Latex');
  t = text(.35,.85,...
           sprintf('$\\propto \\exp(-\\frac{1}{2}(\\frac{X}{%.2f})^2)$', data.fgauss1.s),...
           'Color','k','fontsize',20,'horizontal','left','units','normal');
  set(t,'Interpreter','Latex');
  set(0,'DefaultAxesFontSize','factory');
  set(0,'DefaultTextFontSize','factory');
end


