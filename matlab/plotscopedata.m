clear all 
% close all

% scopefile = 'Axis2_Pulses2.csv';
scopefile = uigetfile('*.csv');

dd = dlmread(scopefile,',',5,0);

plotyy(dd(:,1),dd(:,3),dd(:,1),dd(:,4))
tth = title(scopefile);
set(tth,'interpreter','none')
xlabel('seconds')

legend('CH1','CH2')

L = dd(end,1) - dd(1,1);
Fs = 1/dd(2,1)-dd(1,1);
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(dd(:,3),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')