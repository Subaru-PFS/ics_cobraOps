%Plot the Motormap Correction Curve
%Execute go.m first!
plotMotorMapLineFit(q1.eva.p1(:,:), 41);
plotMotorMapLineFit(q2.eva.p1(:,:), 42);
plotMotorMapLineFit(q3.eva.p1(:,:), 43);
plotMotorMapLineFit(q4.eva.p1(:,:), 44);
plotMotorMapLineFit(q5.eva.p1(:,:), 45);

 limy = [-1 1]
 q1p1m = mean(q1.eva.p1(:,:),1);
   figure(11) 
   
  
   ylim(limy)
   %Throw out outliers:
   fm1  = q1p1m(find(imag(q1p1m) < 1));
   fm1  = fm1(find(imag(fm1) > -1));
    plot(fm1, 'r.');
   coeffs = polyfit(real(fm1),imag(fm1),3);
    
   % Get fitted values
fittedX = linspace(0, 2*pi, 100);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 3);

    
    
   figure(12) 
    plot(mean(q2.eva.p1(:,:),1), 'r.');
    ylim(limy)
  % plot(q2.eva.p1)
   
    figure(13)
    
     plot(mean(q3.eva.p1(:,:),1), 'r.');
   ylim(limy)
%plot(q3.eva.p1)
 
      figure(14)
      plot(mean(q4.eva.p1(:,:),1), 'r.');
        ylim(limy)
 %  plot(q4.eva.p1)
    figure(15) 
    plot(mean(q5.eva.p1(:,:),1), 'r.');
    ylim(limy)
%plot(q5.eva.p1)
 