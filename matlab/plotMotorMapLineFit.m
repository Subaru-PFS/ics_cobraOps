function [ output_args ] = plotMapLineFit( evaMatrix, figNum )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 limy = [-1 1]
 q1p1m = mean(evaMatrix(:,:),1);
   figure(figNum) 
   
  

   %Throw out outliers:
   fm1  = q1p1m(find(imag(q1p1m) < 1));
   fm1  = fm1(find(imag(fm1) > -1));
    plot(fm1, 'r.');
       ylim(limy)
   coeffs = polyfit(real(fm1),imag(fm1),3);
    
   % Get fitted values
fittedX = linspace(0, 2*pi, 100);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
hold on;
plot(fittedX, fittedY, 'r-', 'LineWidth', 3);
saveas(gcf,strcat('motorMapCorrection_',num2str(figNum)),'fig');


end

