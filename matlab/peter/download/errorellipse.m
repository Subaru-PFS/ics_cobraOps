function [x,y]=errorellipse(x1,y1)
%This program calculates necessary elements for plotting an error ellipse

n=length(x1);
%calculate the covariance from standard errors
    covariance=cov(x1,y1)/n; 
% eigenvalues and eigenvectors
    [eigvector,eigvalue]=eig(covariance);
%the length of one of the axes of the "1 sigma error ellipse"
    axis1=sqrt(eigvalue(1,1)); 
    axis2=sqrt(eigvalue(2,2));

% The angle of the ellipse relative to x-axis
    theta=atan(eigvector(2,2)/eigvector(2,1));

% Plot the error ellipse
    p=0:(2*pi)/100:2*pi;
        x=axis2.*cos(p).*cos(theta)-axis1.*sin(p).*sin(theta)+mean(x1);
        y=axis2.*cos(p).*sin(theta)+axis1.*sin(p).*cos(theta)+mean(y1);