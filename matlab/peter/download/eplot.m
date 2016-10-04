%*****************************************************************************************************
%NAME: eplot.m
%AUTHOR: Andri M. Gretarsson
%DATE: 01/21/98
%
%SYNTAX: eplot(X,Y, 'colour')
%
%'colour' defaults to 'b'
%
%This function acts just like the built-in function 'plot' but plots error-bars. The error-bars
%are plotted in the colour given by 'colour'.  'X' and 'Y' are Nx2 matrixes, the first column
%representing the values of the coordinate (x or y), the second column representing the uncertainty
%('error') of those values.  'colour' is a one-letter string which must be one of the letters allowed
%in the built-in matlab function 'plot', to specify the plot colour.  Note that the function exits
%with "hold" set to "off".
%
%This function does not print points in addition to the error bars.  Where the error bars cross, is
%the coordinate point.  This means that if both error bars are exceedingly small complared to the
%coordinate values, the mark will be correspondingly small.  In such situations, it may be better to
%use 'plot' directly and specify that the error is smaller than the size of the mark.
%
%EXAMPLE:
%
%X=[1.0	0.2
%   2.0	0.2]
%Y=[1.0 	0.25
%   2.0	0.25]
%eplot(X,Y,'g')
%
%plots a green cross of width 0.2 and height 0.25 at coordinate (1.0,1.0), and a cross of width 0.2
%and height 0.25 at coordinate( (2.0,2.0).
%
%LAST MODIFIED:  02/16/06
%*****************************************************************************************************

function eplot(x,y,colourstring)

%phm addition 2/16/06
if ~exist('colourstring')
  colourstring = 'b';
end
if (length(colourstring) > 1) 
  disp('EPLOT: third argument must be a single color argument');
  return;
end
if isempty(regexp(colourstring,'[rgbcmykw]'))
  disp('EPLOT: valid color arguments are r,g,b,c,m,y,k,w');
  return;
end
init_nohold = ~ishold;

xvalue=x(:,1); %For clarity
xerror=x(:,2);
yvalue=y(:,1);
yerror=y(:,2);

plot(xvalue-xerror,yvalue-yerror,'w-',xvalue+xerror,yvalue+yerror,'w-');    
hold on;
%Sets appropriate axes but otherwise invisible on a white background

for i=1:length(xvalue)
   plot([xvalue(i)-xerror(i) xvalue(i)+xerror(i)], ...
	[yvalue(i) yvalue(i)], colourstring);
   plot([xvalue(i) xvalue(i)], ...
	[yvalue(i)-yerror(i) yvalue(i)+yerror(i)], colourstring);
end
if init_nohold
  hold off;
end