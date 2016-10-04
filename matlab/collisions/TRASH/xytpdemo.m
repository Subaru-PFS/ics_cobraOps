function xytpdemo(theta0,d_cobra)
% demonstrate x,y to theta, phi transformation for a cobra
% phi is the phi-arm angle wrt the theta arm
% r_err is the radius of the error circle, also in arm-lengths!
% (.0021 == 5 microns)
  psym = 'r';
  arm = 2.375;
  rfib = 1.0; % fiber arm, elbow radius
  rpatrol = 2*arm;
  %d_cobra = 8.0; % cobra spacing
  
  
  phi_angles = linspace(pi - asin((8 - 2 - arm)/(2*arm)),0,30);
  
  for phi = phi_angles

    % elbow(1) and fiber(2) xy positions
    xy(1) = arm * exp(i*theta0) + d_cobra;
    xy(2) = arm * exp(i*(theta0 + phi)) + xy(1);
    side1 = xy + rfib * exp(i*(theta0+phi+pi/2));
    side2 = xy + rfib * exp(i*(theta0+phi-pi/2));
    
    % elbow and fiber circles in xy space
    elbowXY = xy(1) + rfib * exp(i*2*pi*(0:.01:1));
    fiberXY = xy(2) + rfib * exp(i*2*pi*(0:.01:1));
    J2XY = [linspace(side1(1),side1(2),30), linspace(side2(2),side2(1),30)];

    % filter out fiberXY's that are outside rpatrol
    elbowXYin = elbowXY(abs(elbowXY) < rpatrol);
    fiberXYin = fiberXY(abs(fiberXY) < rpatrol);
    J2XYin = J2XY(abs(J2XY) < rpatrol);

    % calculate TP for elbow and fiber
    elbowTP = XY2TP(elbowXYin, arm);
    fiberTP = XY2TP(fiberXYin, arm);
    J2TP = XY2TP(J2XYin, arm);

    subplot(221)
    plot(elbowXY,'b--');
    hold on;
    plot(fiberXY,'b');
    plot(J2XY,'b');
    %  hold off;
    axis equal;
    plotcircle(0,0,2*arm,'k');
    xlabel('X [mm]')
    ylabel('Y [mm]')
    title('PFI bench coordinates (X,Y)')
    
    subplot(222)
    plot(elbowTP.tht/pi,elbowTP.phi/pi,'b--','markersize',6);
    hold on;
    plot(fiberTP.tht/pi,fiberTP.phi/pi,'b','markersize',6);
    plot(J2TP.tht/pi,J2TP.phi/pi,'b');
    axis equal;
    xlabel('\theta/\pi');
    ylabel('\phi/\pi');
    title('Cobra coordinates (\theta,\phi) [detail]');
    
    
    subplot(212)
    plot(fiberTP.tht/pi,fiberTP.phi/pi,'b','markersize',1);
    hold on;
    plot(elbowTP.tht/pi,elbowTP.phi/pi,'b--','markersize',1);
    plot(J2TP.tht/pi,J2TP.phi/pi,'b');
    xlim([-1 1]);
    ylim([0 1]);
% $$$   axis equal;
    xlabel('\theta/\pi');
    ylabel('\phi/\pi');
    title('Cobra coordinates (\theta,\phi)');
    drawnow;

  end
  thick = 3;
  subplot(221)
  plot(elbowXY,'k','linewidth',thick);
  plot(fiberXY,'k','linewidth',thick);
  plot(J2XY,'k','linewidth',thick);
  %  hold off;
  
  subplot(222)
  plot(elbowTP.tht/pi,elbowTP.phi/pi,'k','linewidth',thick);
  hold on;
  plot(fiberTP.tht/pi,fiberTP.phi/pi,'k','linewidth',thick);
  plot(J2TP.tht/pi,J2TP.phi/pi,'k','linewidth',thick);
  
  subplot(212)
  plot(fiberTP.tht/pi,fiberTP.phi/pi,'k','linewidth',thick);
  plot(elbowTP.tht/pi,elbowTP.phi/pi,'k','linewidth',thick);
  plot(J2TP.tht/pi,J2TP.phi/pi,'k','linewidth',thick);
end

function TP = XY2TP(xy, arm)
  R = abs(xy);
  T = angle(xy);
  phi = acos((R/arm).^2*0.5 - 1);
  tht = mod(T - acos(R/arm * 0.5)+pi, 2*pi)-pi;
  TP = tht + i*phi;
end
  
% $$$ 
% $$$ for x0 = z
% $$$   x = x0;
% $$$   y = 0;
% $$$   r = min(r_err, 2-abs(x0));
% $$$   
% $$$   
% $$$   rr = abs(x+i*y);
% $$$   tt = angle(x+i*y);
% $$$   pp = acos(rr^2*0.5 - 1);
% $$$   th = mod(tt - acos(rr/2) - theta0, 2*pi);
% $$$   
% $$$   c_xy = x + i*y + r*exp(i*2*pi*(0:.01:1));
% $$$   R = abs(c_xy);
% $$$   T = angle(c_xy);
% $$$   phi = acos(R.^2 * 0.5 - 1);
% $$$   theta = mod(T - acos(R/2)-theta0, 2*pi); % note that acos(R/2) only
% $$$                                            % works for equal arms
% $$$   
% $$$ % $$$   plot((theta-th),(phi-pp));
% $$$ % $$$   axis equal;
% $$$ % $$$   xlabel('\Delta\theta [rad]');
% $$$ % $$$   ylabel('\Delta\phi [rad]');
% $$$ % $$$   hold on;
% $$$   
% $$$   %  fprintf(1,'%f: %f \n', (pi - pp)*180/pi, x0/(max(theta) - min(theta)));
% $$$   
% $$$   if (length(z) > 100)
% $$$     subplot(221)
% $$$     plot(c_xy*arm,psym);
% $$$     hold on;
% $$$     plot([complex(0), arm*exp(i*(th + theta0)), ...
% $$$           arm*exp(i*(th + theta0)) + arm*exp(i*(th + theta0 + pp))]);
% $$$     hold off;
% $$$     axis equal;
% $$$     plotcircle(0,0,2*arm,'k');
% $$$     xlabel('X [mm]')
% $$$     ylabel('Y [mm]')
% $$$     title('PFI bench coordinates (X,Y)')
% $$$     
% $$$     subplot(222)
% $$$     plot(theta/pi,phi/pi,'r','markersize',6);
% $$$     axis equal;
% $$$     xlabel('\theta/\pi');
% $$$     ylabel('\phi/\pi');
% $$$     title('Cobra coordinates (\theta,\phi) [detail]');
% $$$     
% $$$     
% $$$     subplot(212)
% $$$     plot(theta/pi,phi/pi,'r.','markersize',1);
% $$$     xlim([0 2]);
% $$$     ylim([0 1]);
% $$$ % $$$   axis equal;
% $$$     xlabel('\theta/\pi');
% $$$     ylabel('\phi/\pi');
% $$$     title('Cobra coordinates (\theta,\phi)');
% $$$     drawnow;
% $$$   end
% $$$   %  pause(1)
% $$$ end
% $$$ hold off;
% $$$ 
% $$$ % $$$ if length(z) == 1
% $$$ % $$$   subplot(121)
% $$$ % $$$   plot(c_xy*arm,psym);
% $$$ % $$$   hold on;
% $$$ % $$$   plot([complex(0), arm*exp(i*(th + theta0)), ...
% $$$ % $$$         arm*exp(i*(th + theta0)) + arm*exp(i*(th + theta0 + pp))]);
% $$$ % $$$   hold off;
% $$$ % $$$   axis equal;
% $$$ % $$$   plotcircle(0,0,2*arm,'k');
% $$$ % $$$   xlabel('X [mm]')
% $$$ % $$$   ylabel('Y [mm]')
% $$$ % $$$   title('PFI bench coordinates (X,Y)')
% $$$ % $$$   
% $$$ % $$$   subplot(122)
% $$$ % $$$   plot((theta-th),(phi-pp),'r');
% $$$ % $$$   %  axis equal;
% $$$ % $$$   xlabel('\Delta\theta');
% $$$ % $$$   ylabel('\Delta\phi');
% $$$ % $$$   %  title('Cobra coordinates (\theta,\phi) [detail]');
% $$$ % $$$   
% $$$ % $$$ end