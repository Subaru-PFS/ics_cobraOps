% more concrete examples of how arrays are preferable to cell
% arrays

% run this from Dropbox/PFS_EM/TEST_RESULTS/Metrology

close all
clear all

load phiMtrlgySummary_051414

  disp('%% here are the variables loaded from ''phiMtrlgySummary_051414'':');
  whos

q = input('hit <enter> to continue...');

  disp('%% turn cell array ''data'' into an array');
  disp('>> dmat = cell2mat(data)');
  dmat = cell2mat(data)

q = input('hit <enter> to continue...');

  disp('%% let''s see what the fields are...');
  disp('>> dmat(1)');
  
  dmat(1)
  
q = input('hit <enter> to continue...');

  disp('%% ''cBar'' looks like an array of fitted center positions');
  disp('%% Check out the plot of the cBar positions...');
  
  disp('cbar = vertcat(dmat.cBar);'); cbar = vertcat(dmat.cBar);
  figure(1);
  disp('plot(cbar,''o'')'); plot(cbar,'o');
  
  disp('%% notice the outliers corresponding to cobras 4 and 5')
  
q = input('hit <enter> to continue...');
  disp('%% we can identify them using imagesc');
  
  disp(['>> imagesc(abs(bsxfun(@minus, dmat.cBar, ' ...
        'mean(dmat.cBar))))']);
  figure(2);
  imagesc(abs(bsxfun(@minus, cbar, mean(cbar))))
  
q = input('hit <enter> to continue...');

  disp('%% link length information can also be quickly pulled out');
  
  disp('L2s = vertcat(dmat.Link2s);');
  L2s = vertcat(dmat.Link2s);
  
  figure(3);
  disp('plot(L2s);');
  plot(L2s);