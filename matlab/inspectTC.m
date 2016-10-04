function output = inspectTC(data,cal,minsteps)
% output = inspectTC(data, cal, minsteps);
% data: *_str data structure from Chaz
% cal : structure containing theta0 [deg], center [pix], arm1 [pix],
%       arm2 [pix].
% minsteps: ignore convergences taking fewer steps than minsteps
%           defaults to 1
%
% output: for now, it's an array of the number of steps per convergence.
  
  if ~exist('minsteps','var'), minsteps = 1; end;
  DtoR = pi/180;
  pix2um = 54.37; % ugh! hard coded...
  theta0 = cal.orientation * DtoR;
  cobracenter = cal.s1Center;
  arm1 = cal.L1 * pix2um;
  arm2 = cal.L2 * pix2um;
  RR = @(phi) sqrt(cos(phi) * 2 * arm1 * arm2 + arm1^2 + arm2^2);
  J1err = [data.J1err];
  Rtgt = RR(pi - [data.J2_t]);
  RdJ1 = bsxfun(@times, J1err, Rtgt);
  RdJ2 = arm2 * [data.J2err];
  J2_t = [data.J2_t];
  RdJ1J2 = RdJ1 + i*RdJ2;
  
  for jj=1:length(data)
    e5  = CobraTPerror( 5,pi - data(jj).J2_t,arm1,arm2);
    e10 = CobraTPerror(10,pi - data(jj).J2_t,arm1,arm2);

    in5  = inpoly( e5, RdJ1J2(:,jj));
    in10 = inpoly(e10, RdJ1J2(:,jj));
    try
      nsteps(jj) = find(in5,1);
    catch
      nsteps(jj) = 16;
    end;
    if nsteps(jj) > minsteps
      
      subplot(221)
      plot([0, arm1*exp(i*(data(jj).J1_t + theta0)), ...
            arm1*exp(i*(data(jj).J1_t + theta0)) ...
            + arm2*exp(i*((data(jj).J1_t - pi + data(jj).J2_t) ...
                          + theta0))], 'o-'); 
      hold on;
% $$$       plot(arm1*exp(i*(data(jj).J1 + theta0)) + ...
% $$$            arm2*exp(i*((data(jj).J1 - pi + data(jj).J2) + theta0)),'k.:')
      plot((data(jj).curPos(~isnan(data(jj).status))    - cobracenter)*pix2um,'k.--');
      plot((data(jj).targCmplx - cobracenter)*pix2um,'r.');
      plotcircle(0,0,arm1+arm2,'r');
      %      cmplx(@plotcircle,(data(jj).targCmplx - cobracenter)*pix2um,10,'r');
      cmplx(@plotcircle,(data(jj).targCmplx - cobracenter)*pix2um, 5,'m');
      xlabel('X [um]');
      ylabel('Y [um]');
      title(sprintf('%d: XY coordinates',jj));
      hold off;
      
      subplot(223)
      plot(RdJ1J2(3:end,jj),'ko--');axis equal; hold on;
      cmplx(@scatter,RdJ1J2(3:end,jj),30,(3:15)','fill');
      try, caxis([3 find(data(jj).status,1)]); end;
      colormap('gray');
      cmplx(@text,RdJ1J2(3:end,jj)+1,num2cell(3:15));
      %plot(e10,'r');
      plot(e5,'m');
      xlabel('R\Delta\theta [um]');
      ylabel('r_2\Delta\phi [um]');
      title(sprintf('%d in \\Delta\\theta-\\Delta\\phi space',jj));
      hold off;

      subplot(222)
      semilogy(abs(data(jj).J1err),'o-');hold on;
      plot(abs(data(jj).J2err),'ro-');
% $$$     refline(0,4.2e-3,0,'k--');
      try, refline(find(in5,1),1,Inf,'r-'); end
      %try, refline(find(in10,1),1,Inf,'k--'); end
      %try, refline(find(data(jj).status,1),1,Inf,'g--'); end
      hold off;
      ylabel('angular error [rad]');
      xlabel('step #');
      title('angular error vs. step #')
      legend('\theta','\phi');
      
      subplot(224)
      semilogy(abs(RdJ1(:,jj)),'o-');hold on;
      plot(abs(RdJ2(:,jj)),'ro-');
% $$$     refline(0,10,0,'k--');
      try, refline(find(in5,1),1,Inf,'r-'); end
      %try, refline(find(in10,1),1,Inf,'k--'); end
      %try, refline(find(data(jj).status,1),1,Inf,'g-'); end
      hold off;
      ylabel('distance error [\mum]');
      xlabel('step #');
      title('distance error vs. step #')
      drawnow;
% $$$       pause(2);
      input('next..');
    end
  end
  output = nsteps;
  
end

function output=CobraTPerror(maxerr, J2_t, arm1, arm2)
% calculate a locus of points in dtheta-dphi space corresponding to
% an error circle of radius maxerr in XY space
%
% maxerr: radius of error circle in XY space in MICRONS
% J2_t: target joint 2 angle in RADIANS
% arm1: J1 arm length in MICRONS, defaults to 2375
% arm2: J2 arm length in MICRONS, defaults to 2375

  if ~exist('arm1','var'), arm1 = 2375; end
  if ~exist('arm2','var'), arm2 = 2375; end
  
  if J2_t > pi
    warning('J2_t should be in radians, not degrees');
  end
  
  RR = sqrt(cos(J2_t) * 2 * arm1 * arm2 + arm1^2 + arm2^2);
  
  % for simplicity, assume target lies on the Y axis (X = 0)
  xyerror = i*RR + maxerr * exp(i*2*pi*(0:.01:1));

  J1_t = pi/2 - asin(arm2 * sin(J2_t) / RR);
  
  % force error bounds to be within the patrol region
  xyR = min(abs(xyerror), arm1+arm2);
  xyT = angle(xyerror);
  xyerror = xyR .* exp(i*xyT);

  PHI = acos((xyR.^2 - arm1^2 -arm2^2) / (2 * arm1 * arm2));
  THT = xyT - asin(arm2 * sin(PHI) ./ xyR);
  
  output = RR*(THT - J1_t) + i*arm2*(PHI - J2_t);
  output = output(:);
end