function [output  pltHndl] = calcSNR(data,badthresh,serialNum,t_step,t_max,t_obs,k_offset)
% [output pltHndl] = calcSNR(data, badthreshold, [serialNum, t_step, t_max, t_obs, k_offset]);
% data: *_str data structure from Chaz
% serialNum: serial number or -1 for auto-pId (-1)
% t_step: time per step [s] (8)
% t_max: maximum configuration time [s] (105)
% t_obs: baseline observation time [s] (900)
% k_offset: loss parameter [1/mm^2] (1/(.075)^2): loss = k dx^2
%
% output: 
%  output.snr        : (#iterations,#targets) array of SNR's (min value 0)
%  output.maxsnr     : maximum snr averaged over targets
%  output.indx_maxsnr: iteration index of maxSNR (matlab indexing)
%  output.bad        : target indices of poor performers


% by: Peter Mao

    if ~exist('badthresh','var'), badthresh = 0.9; end;
    if ~exist('serialNum','var'),serialNum =   -1; end;  
    if ~exist('t_max','var'), t_max  = 105; end;
    if ~exist('t_obs','var'), t_obs  = 900; end;
    if ~exist('k_offset','var'), k_offset = 0.075^(-2); end;
    DtoR = pi/180;
    pID = data.pid;
    
    firstnan = find(isnan(data(jj).dist), 1);
    dist = data(jj).dist * 1e-3;
  if ~exist('t_step','var')
    t_step_old =  [20,28,36,44,52,60,68,76,84,92,100]';
    t_step = [12,20,28,36,44,52,60,68,76,84,92]';
  end;
    if length(dist)==1 || isempty(dist)
        continue;
    end
    dist(firstnan:end) = dist(firstnan - 1);
 try
     timeOld =  t_max + t_obs - t_step_old;
    time =  t_max + t_obs - t_step; 
    output(jjc).snr = (1 - dist.^2 * k_offset) .* sqrt(time / t_obs);
    output(jjc).snrOld = (1 - dist.^2 * k_offset) .* sqrt(timeOld / t_obs);

     catch
     keyboard
 end
    negs = output(jjc).snr < 0;
    output(jjc).snr(negs) = 0;
    jjc=jjc+1;
  end

%   keyboard;
  [maxsnr indx_maxsnr] = max(mean([output.snr],2));
  t_config = t_obs+t_max - time;
%   output.maxsnr = maxsnr;
%   output.indx_maxsnr = indx_maxsnr;

if serialNum == -1
    titlestring = sprintf('pID %d, snr= %.3f at %d',pID,maxsnr,indx_maxsnr-1);
else
    titlestring = sprintf('EM-%d, snr= %.3f at %d',serialNum,maxsnr,indx_maxsnr-1);
end

  %%%%%%summary plot
  plot(t_config, [output.snr],'.'); hold on;
  pltHndl = plot(t_config, mean([output.snr],2), 'o-', 'Linewidth',4,'MarkerSize',5);
%   text(t_config+1, mean([output.snr],2)+.01,...
%        {'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14'},...
%        'fontsize',14)
%   text(t_config+1, mean([output.snr],2)+.01,...
%        {'0','1','2','3','4','5','6','7','8','9'},...
%        'fontsize',14)
  hold off;
  ylim([.8 1.05]);
  grid on;
  xlabel('configuration time');
  ylabel('SNR/SNR_{900s}');
  title(titlestring,'fontsize',12);

  %% text output
  fprintf(1,'max SNR = %.3f at %d iterations\n',maxsnr, indx_maxsnr - 1);