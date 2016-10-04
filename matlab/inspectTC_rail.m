function output=inspectTC_rail(field)
% usage: inspectTC_rail(field_number)
%
% run this in the directory of interest.
%
% generates visualizations of convergence of cobras over the whole
% rail.

%CONSTANTS
  if ~exist('field','var'), field = 1; end;
  um_per_pixel = 90.0;
  R_tip = 1000 / um_per_pixel; % tip radius in pixels.

  data0 = loadmats;
  
  % move all the data structures with names ending in _str into a
  % numbered cell array for easy access
  kk = 0;
  for ff=sort(fields(data0)')
      if ~isempty(regexp(ff{1},'_str$'))
         kk = kk+1;
         pdata{kk} = data0.(ff{1});
         pids(kk) = pdata{kk}(1).pid;
      end
  end
  nfields = length(pdata{kk});
  
  
  %% read in the calibration file
  geom = defineBenchGeometry([],1,1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% overwrite geom structure entries to remove fibers that are not relevant to this data set
  [tmp, gKeep] = ismember(pids,geom.pids);
  [nnKeepR nnRow] = ismember(geom.NN.row,gKeep);
  [nnKeepC nnCol] = ismember(geom.NN.col,gKeep);
  nnKeep = nnKeepR & nnKeepC;
  geom.NN.row = nnRow(nnKeep);
  geom.NN.col = nnCol(nnKeep);
  geom.NN.xy  = geom.NN.xy(nnKeep);
  
  nGeomPids = length(geom.pids);
  for ff=fields(geom)'
      fld = ff{1};
      if (size(geom.(fld),1) == nGeomPids)
          if (size(geom.(fld),2) == nGeomPids)
              % #pid X #pid square matrices
              geom.(fld) = geom.(fld)(gKeep,gKeep);
          else % #pid X anything non-square matrices
              geom.(fld) = geom.(fld)(gKeep,:);
          end
      end
  end
  
  geom.alpha = 0; % errors are already baked in, we don't need to
                  % add any.
                  %  geom.thteps = .0001;
  clear data0 tmp nnKeepR nnKeepC

  %% GEOM READY FOR USE: index in geom agrees with data index
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  unpack some data from geom
  center = geom.center.';
  L1 = geom.L1.';
  L2 = geom.L2.';
  tht0 = geom.tht0.';
  phiIn = geom.phiIn.';
  
  %% plot control junk
  curaxlim = [];
  %%
  
  disp(['Enter (without quotes): ''-'' to decrement, '''' to increment, ''#'' to go directly ' ...
        'to target #, ''nan'' to finish'])
  
  while ~isnan(field)

      % pull the data for the requested field (all fibers)
      for kk=1:length(pdata)
          J1(:,kk)   = pdata{kk}(field).J1 + geom.tht0(kk);
          J2(:,kk)   = pdata{kk}(field).J2 - pi;
          tgt(kk)    = pdata{kk}(field).targCmplx;
          try
          j1t     = pdata{kk}(field).J1_t(5) + geom.tht0(kk);
          j2t     = pdata{kk}(field).J2_t(5) - pi;
          catch % backward compatibility before targets were a list. 
          j1t     = pdata{kk}(field).J1_t + geom.tht0(kk);
          j2t     = pdata{kk}(field).J2_t - pi;
          end
          tgt2(kk)   = L1(kk) * exp(i*j1t) + L2(kk) * exp(i*(j1t+j2t)) + center(kk);
      end
      
      disp(mod(angle(tgt - center) - angle(tgt2 - center), 2*pi)/pi)
      % tip positions
      xy = (bsxfun(@plus,center,bsxfun(@times, L1, exp(i*J1))) + ...
            bsxfun(@times, L2, exp(i*(J1+J2)))).';
      
      %% generate the plot
      clf;
      % patrol area
      cmplx(@plotcircle,center,L1+L2,'k');   hold on;
      cmplx(@text,center,cellstr(num2str(pids(:)))); % PID labels for regions
                                                     % tht0 direction
      plot([center; center + (L1+L2) .* exp(i*tht0)],'k--');
      % home position
      plot(center + L1 .* exp(i*tht0) + L2 .* exp(i*(tht0 + phiIn)), 'go','MarkerFace','g'); 
% $$$   % tip @ home (full extent)
% $$$   cmplx(@plotcircle, center + L1 .* exp(i*tht0) + L2 .* exp(i*(tht0 + phiIn)), R_tip, 'g');
% $$$   % elbow @ home (full extent)
% $$$   cmplx(@plotcircle, center + L1 .* exp(i*tht0), R_tip,'k'); 
% target location
      plot(tgt,'ro','Markerface','r');
      plot(tgt2,'rx','MarkerSize',10); 
      set(gca,'xdir','reverse')


% $$$   disp('using perfect moves');
% $$$   geom.S1Pm = [];
      
      xystart = xy(:,1);
      for tt=1:9 % loop over time
          trj = [];
          thisTraj = generateTrajectory(xystart,xy(:,tt+1),geom,'earlyLate');
          trj = thisTraj.traj;
          xystart = trj(:,end); 
          trjTP = XY2TP(bsxfun(@minus,trj,geom.center), geom.L1, geom.L2);
          elbow = bsxfun(@plus,geom.center,...
                         bsxfun(@times, geom.L1, exp(i*trjTP.tht)) );
          plot(trj.'); % fiber trajectory
          plot(elbow.','--'); % elbow trajectory
          
          plot([geom.center elbow(:,1) trj(:,1)].','b')
          plot([geom.center elbow(:,end) trj(:,end)].','b') 
      end
      title(sprintf('target #%d',field));
      if ~isempty(curaxlim), axis(curaxlim); end;
      hold off;
      
      num_input = input(sprintf('next target #[%d]: ',mod(field,nfields)+1),'s');
      if isempty(num_input)
          field = field + 1;
      elseif num_input == '-'
          field = field - 1;
      else
          field = str2num(num_input);
      end
      field = mod(field-1,nfields)+1;
      curaxlim = axis;
  end
  
  output=packstruct(geom,J1,J2,xy,tgt);