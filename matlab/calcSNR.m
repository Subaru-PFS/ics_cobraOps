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

    %%%%%% SNR vs time summary plot
    figure(747)
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
    
    %%% SNR vs target
    figure(757)
    plot(0:length(snr)-1, snr, 'b.'); hold on;
    plot(0:length(snr)-1, snr(indx_maxsnr,:), 'ro'); hold off;
    refline(0,badthresh,0,'k:');
    xlabel('target ID (from 0)');
    ylabel('SNR');
    grid on;
    title(titlestring,'fontsize',12);
    
    
    %% text output
    fprintf(1,'max SNR = %.3f at %d iterations\n',maxsnr, indx_maxsnr - 1);

    
    %% run inspectTC_CIT on bad ones
    if serialNum == -1 %%& strcmp(getenv('USER'),'petermao')
% $$$         inspectTC_CIT(pID,1,output.bad)

% $$$         j1  = [data.J1];   % theta
% $$$         j1s = [data.J1_S]; % theta steps commanded
% $$$         dj1 = [data.J1_t] - j1; % angular displacement to target (should
% $$$                                 % correlate with j1s)
% $$$         figure(990)
% $$$         histogram(mod(j1(1,:)+pi,2*pi)-pi,-1:.025:1);
% $$$         title('Initial \theta distribution over [-pi,pi)');
% $$$         for step = 1:1
% $$$             figure(888)
% $$$             plot(j1(step,:),j1s(step,:),'o'); hold on;
% $$$             plot(j1(step,output.bad),j1s(step,output.bad),'ro');hold off;
% $$$             ylabel('steps');
% $$$             xlabel('\theta [rad]');
% $$$             title(sprintf('pre-move %d steps vs. initial angle',step));
% $$$             refline(2*pi,0,Inf,'k:');
% $$$             grid on;
% $$$             figure(889)
% $$$             plot(dj1(step,:), j1s(step,:),'o'); hold on;
% $$$             plot(dj1(step,output.bad), j1s(step,output.bad),'ro'); hold off;
% $$$             set(gca,'XAxisLocation','origin');
% $$$             set(gca,'YAxisLocation','origin');
% $$$             ylabel('steps');
% $$$             xlabel('\Delta \theta [rad]');
% $$$             title(sprintf('pre-move %d steps vs. \\Delta angle',step));
% $$$ % $$$             keyboard;
% $$$         end
    end