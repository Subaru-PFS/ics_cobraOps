function [output posmax jjmax] = calcRailSNR(t_step,t_max,t_obs,k_offset)
% output = calcSNR(data, cal, [t_step, t_max, t_obs, k_offset]);
% 
% t_step: time per step [s] (8)
% t_max: maximum configuration time [s] (105)
% t_obs: baseline observation time [s] (900)
% k_offset: loss parameter [1/mm^2] (1/(.075)^2): loss = k dx^2
%
% output: n-positioners X m-fields snr array of struct [obsolete]
% output: positioner x target x iteration
% posmax: max mean-over-fields snr for position
% jjmax:  iteration index of max mean-over-for
%  
% blame: Peter Mao
  
  if ~exist('t_step','var'),t_step = [12,20,28,36,44,52,60,68,76,84,92]'; end;
  if ~exist('t_max','var'), t_max  = 105; end;
  if ~exist('t_obs','var'), t_obs  = 900; end;
  if ~exist('k_offset','var'), k_offset = 0.075^(-2); end;
  DtoR = pi/180;

  data = []; % _str data structures in a mxn array)
  mats = loadmats;
  for f = sort(fields(mats))'
      data = vertcat(data, mats.(f{1}));
  end

  pID = data.pid;

  datasize = size(data);

  %%%%%% calculation of bkg dominated snr.
  for kk = 1:datasize(2) % iterations, columns
    for jj=1:datasize(1) % positioners, rows
      
      firstnan = find(isnan(data(jj,kk).dist), 1);
      dist = data(jj,kk).dist * 1e-3; % convert to mm
      %% invalid dist data shows up as -1, so increase by 1e6X
      %% for negative dists
      dist(dist < 0) = -1e6;     
      
      dist(firstnan:end) = dist(firstnan - 1);

      time =  t_max + t_obs -  t_step;
% $$$       output(jj,kk).snr = (1 - dist.^2 * k_offset) .* sqrt(time / t_obs);
      output(jj,kk,:) = ((1 - dist.^2 * k_offset) .* sqrt(time / t_obs)).';
% $$$       negs = output(jj,kk).snr < 0;
% $$$       output(jj,kk).snr(negs) = 0;
    end
  end
  output(output < 0) = 0;

  % calculate for each positioner over all fields
  for jj = 1:datasize(1)
% $$$     thispositioners_meanSNR = mean([output(jj,:).snr],2); % at each step
      thispositioners_meanSNR = mean(squeeze(output(jj,:,:)),2); % at each step
      [posmax(jj) jjmax(jj)] = max(thispositioners_meanSNR);
  end
  
  return;
  
  [maxsnr indx_maxsnr] = max(mean([output.snr],2));
  t_config = t_obs+t_max - time;

