function scansdf(data)
% SCANSDF: interactively scan a microXAM image
% USAGE: scansdf(data)
%        data is either a data structure from readsdf
%        OR
%        a file with SDF formatted data
%
% run-time inputs: n,p for next or previous row
%                  N,P jumps by 10 rows
%                  no arg runs last case
%                  q to quit
%                  anything else is interpreted as a matlab command

  if ischar(data)
    data = readsdf(data);
  end
  
  nscan = data.header.NumProfiles;
  
  this_row = floor(nscan/2);
  reply = ' ';
  step = 1;
  while ~strcmp(reply,'q')
    y = data.m(this_row,:);
    plot(data.x,y);
    title(['row: ' num2str(this_row)])
    reply = input('N,n,P,p,q> ','s');
    switch reply
     case 'N'
      step = 10;
     case 'n'
      step = 1;
     case 'P'
      step = -10;
     case 'p'
      step = -1;
     case {'','q'}
     otherwise
      try,
        eval(reply),
      catch,
        disp(['command ' reply ' not recognized']),
      end,
      step = 0;
    end
    this_row = this_row + step;
    this_row = min(nscan-2,this_row);
    this_row = max(1,this_row);
  end
      