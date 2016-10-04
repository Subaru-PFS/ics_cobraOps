function output=analyzeTC_Joe(data,cal,maxsteps)
% output=analyzeTC(data,cal,maxsteps)
% data: *_str data structure from Chaz
% cal : pId# structure from ../../../SVN/MATLAB/EMconfig_MM-DD-YY.mat
% maxsteps: ignore convergences taking more steps than maxsteps
%
% output: J1, J2 error/request arrays,
% if REDUCE works, J1r and J2r hold the alpha and beta values for
% the motors.
%
% if the motor is underdamped, there is commented out code to remove
% the damping with a line fit to error vs. request.
  
  if ~exist('maxsteps','var'), maxsteps = 16; end
  
  %% motor calibration constants
  DtoR = pi/180;
  PDE = 2.0; %um
  pid = data(1).pid;
  pix2um = getARMval(cal,pid,'Pixel_scale');
  theta0 = getARMval(cal,pid,'Global_base_ori_z') * DtoR;
  cobracenter = getARMval(cal,pid,'Global_base_pos_x') + ...
      getARMval(cal,pid,'Global_base_pos_y') * i;
  arm1 = getARMval(cal,pid,'Link1_Link_Length') * pix2um;
  arm2 = getARMval(cal,pid,'Link2_Link_Length') * pix2um;
  RR = @(phi) sqrt(cos(phi) * 2 * arm1 * arm2 + arm1^2 + arm2^2);
 

  ntargets = length(data);
  kk = 1;
  %% init outputs to column vectors
  J1.pos = [0;0];
  J1.req = [0;0];  % all moves request
  J1.err = [0;0];  % all moves error
% $$$   J1.req1 = [0;0]; % first move request
% $$$   J1.err1 = [0;0]; % first move error
  J1.iter = [0;0]; % iteration index
  J1.req_pde = [0;0]; % angular error
  J1.err_pde = [0;0]; % angular error
  J1.vec = [0 0 ; 0 0 ];
  J2 = J1;
  BETA = 1; % 1 for johannes, .5 for peter
  
  motorMap = linspace(3.6,360,100);
  
  bucketList = [];
  fwdvalues = [];
  rvsvalues = [];
 
  for jj=1:ntargets % put all datapoints in one big list
      diffv = diff(data(jj).J1);
      for kk=1:length(data(jj).J1_S)-1
         % if(abs(data(jj).J1(kk))<0.01) % filter out first move
          yval = 1i* diffv(kk)/data(jj).J1_S(kk)*180/pi; % Motor map value in rad/step
          val1 = data(jj).J1(kk) + yval;
          val2 = data(jj).J1(kk+1) + yval;
          if(data(jj).J1_S(kk)>0) % Sort fwd and rvs moves.
              %values = vertcat(values, [val1, val2]);
              fwdvalues = vertcat(fwdvalues, [val1, val2]);
          else
              rvsvalues = vertcat(rvsvalues, [val1, val2]);
          end
        %  end
      end
  end
  keyboard;
  %filter out NaNs and Infinite MM Values (last moves) 
   fwdvalues(any(isnan(fwdvalues),2),:)=[];
   fwdvalues = fwdvalues(isfinite(fwdvalues(:, 1)), :);
   fwdvalues(imag(fwdvalues(:,1))<0,:)=[];
   fwdvalues(real(fwdvalues(:,1))<0.02,:)=[];
     %filter out NaNs and Infinite MM Values (last moves) 
   rvsvalues(any(isnan(rvsvalues),2),:)=[];
   rvsvalues = rvsvalues(isfinite(rvsvalues(:, 1)), :);
   rvsvalues(imag(rvsvalues(:,1))<0,:)=[];
  
   figure(11) 
   plot(real(fwdvalues(:,2)),imag(fwdvalues(:,2)),'*r')
   hold on
   plot(real(fwdvalues(:,1)),imag(fwdvalues(:,1)),'*b')
   plot(fwdvalues.');
   hold off
   
   
   figure(12)
   plot(real(rvsvalues(:,2)),imag(rvsvalues(:,2)),'*r')
   hold on
   plot(real(rvsvalues(:,1)),imag(rvsvalues(:,1)),'*b')
   plot(rvsvalues.');
   hold off
   keyboard;
%plot([real(fwdvalues(:,1)); imag(fwdvalues(:,1))*180/pi ],[real(fwdvalues(:,2)); imag(fwdvalues(:,2))*180/pi],'-k')


  for jj=1:ntargets
      stat0 = [find(data(jj).status == 0) ;...
               find(data(jj).status == 1, 1)];
      kkstrt = stat0(1:end-1);
      kkfnsh = stat0(2:end);
      
      nsteps = length(kkstrt);
      kkrange = kk:(kk+nsteps-1);
      Rpos   = RR(data(jj).J2(stat0)); % fiber radial position in
                                       % patrol region
      %% throw the baby out with the bath water (anything with
      %% error > 1.5 radian on either motor gets ejected)
      if (max(abs([data(jj).J1err(kkfnsh) ; data(jj).J2err(kkfnsh)])) < 2) 
        
      % put the data in the buckets
      
      for ii = 1:100 % for all motor map buckets
         
      end
      
    end
  end
  
  % calculate the fraction of the move that is within the bucket
  
  % calculate the expected mm value from existing mm
  
  % multiply the fractions of the move time the mm = 
  
  % weigh move by offset to existing value
  
 
    clear output;
    output = packstruct(fwdvalues, rvsvalues);

  
end
