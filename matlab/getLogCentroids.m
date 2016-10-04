function output = getLogCentroids(msimLogFile)
% msimLogFile is an output_log.txt file
% by default, works in the current directory.
  
  if ~exist('msimLogFile','var'), msimLogFile = 'output_log.txt'; end;
  if ~exist(msimLogFile,'file')
    fprintf(1,'File %s does not exist...exiting',msimLogFile);
    output = [];
    return;
  end;
  logF = fopen(msimLogFile);


  % Line filters
  filedate = regexp(msimLogFile,'\d*_\d*_\d*_\d*_\d*_\d*','match');
  fltSvImg = strcat('seq_saveCurrent_RemoteImage.*',filedate);
  fltTemp = 'Get_Temp';

  % Initiate vars
  tline  = ' ';
  refidx = 1;
  imgidx = 0; % this gets incremented after h87idx.
  sfidx  = 0;
  ffidx  = 0;
  
  RR = exp(-i*(pi/4 - 57/2500));
  
  while ischar(tline)
    % Get next line
    tline = fgetl(logF);
    clear idxHorn;
    idxHorn = strfind(tline,'Horn');
    if ~isempty(idxHorn);
      substrings = strsplit(tline,' ==>> ');
      if length(substrings) == 2
        cmd = substrings{1};
        ret = substrings{2};
        numbers = regexp(ret,'[\d.]+','match');
        clear substrings;
        switch cmd
          case 'Set_HornMethod_FiducialCoordinates()'
            % do nothing.  this is handled by next case anyways.
          case 'cmd_setHornMethodFiducialCoordinate()'
            refpos(refidx) = str2num(numbers{1}) + str2num(numbers{2})*i;
            refidx = refidx + 1;
          case 'Apply_HornMethod()'
            substrings = strsplit(ret, ' ');
            switch substrings{1}
              case 'outRot[3][3]:' % look for data on following line
                if imgidx > 0 % sort sf locations by RailX
                  [tmp idx] = sort(real(sfraw(imgidx,:) * RR));
                  sfraw(imgidx,:) = sfraw(imgidx,idx);
                  sfl0(imgidx,:)  = sfl0(imgidx,idx);
                end
                imgidx = imgidx + 1;
                sfidx = 1;
                ffidx = 1;
                tline = fgetl(logF);
                numbers = regexp(tline,'[\d.]+','match');
                h87(imgidx).R = angle(str2num(numbers{1}) - str2num(numbers{2})*i);
              case 'outT[3]:' % look for data on following line
                tline = fgetl(logF);
                numbers = regexp(tline,'[\d.]+','match');
                h87(imgidx).T = str2num(numbers{1}) + str2num(numbers{2})*i;
              case 'outScale'
                h87(imgidx).S = str2num(numbers{1});
              case 'Weighted'
                rawX = str2num(numbers{2});
                rawY = str2num(numbers{3});
                RailY = imag((rawX + i*rawY) * RR);
                if (myfrac(rawX) ~= 0 | myfrac(rawY) ~= 0)
                  if RailY > -148  % SF's are between -20 and -150 in RailY
                    if RailY < -20
                      sfraw(imgidx,sfidx) = rawX + rawY * i;
                      sfl0(imgidx,sfidx)  = str2num(numbers{4}) + str2num(numbers{5})*i;
                      sfidx = sfidx + 1;
                    elseif abs(RailY - 210) < 5 % FF's are at 210 +-5 in RailY
                      ffraw(imgidx,ffidx) = rawX + rawY * i;
                      ffl0(imgidx,ffidx)  = str2num(numbers{4}) + str2num(numbers{5})*i;
                      ffidx = ffidx + 1;
                    end
                  end
                end
              case 'Fiducial'
                ffraw(imgidx,ffidx) = str2num(numbers{2}) + str2num(numbers{3})*i;
                ffl0(imgidx,ffidx)  = str2num(numbers{4}) + str2num(numbers{5})*i;
                ffidx = ffidx + 1;
              otherwise
                disp(ret);
            end
          otherwise
            % print out oddballs to screen
            fprintf(1,'%s: %s\n',cmd, ret);
        end % switch
      end
    end
  end % while
      % sort sf locations by RailX
  [tmp idx] = sort(real(sfraw(imgidx,:) * RR));
  sfraw(imgidx,:) = sfraw(imgidx,idx);
  sfl0(imgidx,:)  = sfl0(imgidx,idx);

  fclose(logF);

  output.ref    = refpos;
  output.horn   = h87;
  output.sf.raw = sfraw;
  output.sf.L0  = sfl0;
  output.ff.raw = ffraw;
  output.ff.L0  = ffl0;
end

function retval = myfrac(input)
  retval = input - fix(input);
end