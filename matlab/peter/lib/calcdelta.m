function [delta, stdev_delta] = calcdelta(minor_isotope, major_isotope, ...
					 standard, std_minor, std_major )
% CALCDELTA calculates the delta deviation from a given standard
% USAGE: [result stdev_delta] = calcdelta(minor_isotope, major_isotope,standard,
%                           [minor_isotope_stdev], [major_isotop_stdev])
%
% examples: delta17 = calcdelta(O17, O16, smow17)
%           [d18 sigma_d18] = calcdelta(O18, O16, smow18, sig_O18, sigO16)
% if the isotope data are arrays, then result will also be an array

delta = ((minor_isotope./major_isotope)/standard - 1)*1000;

if ( exist('std_minor') && exist('std_major') )
  % the partial derivatives of delta wrt each isotope
  d_minor = major_isotope.^(-1) * (1000/standard);
  d_major = (minor_isotope ./ major_isotope.^2) * (1000/standard);
  stdev_delta = sqrt( (std_minor .* d_minor).^2 + ...
		      (std_major .* d_major).^2 );
else
  stdev_delta = NaN;
  if (nargout == 2)
    disp('standard deviations of the major and minor species must');
    disp('be specified for stdev_delta to be calculated');
  end
end