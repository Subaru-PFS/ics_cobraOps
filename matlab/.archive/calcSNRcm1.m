function [output maxsnr indx_maxsnr] = calcSNR(data,t_step,t_max,t_obs,k_offset)
% output = calcSNR(data, cal, [t_step, t_max, t_obs, k_offset]);
% data: *_str data structure from Chaz
% t_step: time per step [s] (8)
% t_max: maximum configuration time [s] (105)
% t_obs: baseline observation time [s] (900)
% k_offset: loss parameter [1/mm^2] (1/(.075)^2): loss = k dx^2
%
% output: 
%  
% blame: Peter Mao
  
  if ~exist('t_step','var'),t_step =   8; end;
  if ~exist('t_max','var'), t_max  = 105; end;
  if ~exist('t_obs','var'), t_obs  = 900; end;
  if ~exist('k_offset','var'), k_offset = 0.075^(-2); end;
  DtoR = pi/180;
  pID = data.pid;

  %%%%%% calculation of bkg dominated snr.
  for jj=1:length(data)
    
    firstnan = find(isnan(data(jj).dist), 1);
    dist = data(jj).dist * 1e-3;
    dist(firstnan:end) = dist(firstnan - 1);

    time =  t_max + t_obs - (1:length(dist)).' * t_step;
    output(jj).snr = (1 - dist.^2 * k_offset) .* sqrt(time / t_obs);
    negs = output(jj).snr < 0;
    output(jj).snr(negs) = 0;
  end

  [maxsnr indx_maxsnr] = max(mean(vertcat(output.snr),2));
  t_config = t_obs+t_max - time;
%   output.maxsnr = maxsnr;
%   output.indx_maxsnr = indx_maxsnr;

  titlestring = sprintf('pID %d, max(<snr>) = %.3f at %d iterations',pID,maxsnr,indx_maxsnr-1);
  disp(titlestring);

  %%%%%%summary plot
  CW = jet(length(output));
  
  %Create text array of numbers
  maxIters = 15;
  for ii=1:maxIters
      txtArr{ii} = num2str(ii);
  end
  
  for ii=1:length(output)
      plot(t_config, output(ii).snr,'.','color',CW(ii,:)); hold on;
      plot(t_config, mean(output(ii).snr,2), 'o-', 'Linewidth',4,'MarkerSize',5,'color',CW(ii,:));
    %   text(t_config+1, mean(output(ii).snr,2)+.01,...
    %        {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14'},...
    %        'fontsize',14)
      text(t_config(ii)+1, mean(output(ii).snr,2)+.01,...
           txtArr{ii},...
           'fontsize',14)
      hold off;
  end
  ylim([.8 1.05]);
  grid on;
  xlabel('configuration time');
  ylabel('SNR/SNR_{900s}');
  title(titlestring,'fontsize',14);
