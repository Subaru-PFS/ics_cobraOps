function dataout = readcdf(filename,dt_unit)
% MegaSIMS cdf reader
%
% USAGE: data = readcdf(''file[.cdf]'',dt_unit)
%        dt_unit = '[s]ec','[m]in','[h]our','[d]ay'
% 
% defaults: dt_unit = sec
%
% OUTPUT:structure containing cdf data

% --> 1.3 change log
%     make dt.priority=time.priority
% --> 1.2 change log
%     remove multiple epoch support
%     change "rawheader" fieldname to "header"
%     remove hardcoded priority flag
%     add background flag support
%     global variables of type char are left in cell arrays

  ISOCTAVE = logical(exist('OCTAVE_VERSION'));
  FieldSeparationString = '__';
  
  if (~exist('dt_unit'))
    TimeConversionFactor = 86400;
    dt_unit = 's';
  end
  switch dt_unit
   case {'s','sec'}
    TimeConversionFactor = 86400;
   case {'m','min'}
    TimeConversionFactor = 1440;
   case {'h','hour'}
    TimeConversionFactor = 24;
   case {'d','day'}
    TimeConversionFactor = 1;
   otherwise
    disp('Abort [READCDF.M]: illegal time unit for dt. use sec,min or hour')
    dataout = NaN;
    return
  end
  if ISOCTAVE
    TimeConversionFactor = TimeConversionFactor / 86400 / 1000;
  end
  %append '.cdf' suffix if necessary
  if ( (    ISOCTAVE && ~index(filename,'.cdf')    ) ...
       || (~ISOCTAVE && isempty(strfind(filename,'.cdf'))) )
    filename = [filename '.cdf'];
  end
  % THIS IS THE MATLAB SUPPLIED ROUTINE FOR CDF
  if ISOCTAVE
    [data_in header] = ocdfread(filename);
  else
    [data_in header] = cdfread(filename);
  end
  dataout.header = header;
  dataout.header.Path = pwd;
  
  % GET SIZE OF DATA ARRAY, NUMBER OF EPOCH DATA COLUMNS
  [n_datapoints n_variables] = size(data_in);
  if (n_datapoints == 0)
    warning('READCDF:nodata','READCDF: No data in file');
    return;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%% REAL WORK STARTS HERE %%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %extract the FLAG_bkg field
  % determine the location and index of the FLAG_bkg data vector
  for jj=1:n_variables
    if ~isempty(findstr(lower(header.Variables{jj}),'bkg'))
      FLAG_bkg = logical([data_in{:,jj}]);
      COL_bkg = jj;
      if ( length(FLAG_bkg(FLAG_bkg == 0)) == 0 )
	%ie, if there are no data entries
	disp(['Warning [READCDF.M]: No data points found.'])
%	return;
      end
      break;
    end
  end
  if ~exist('COL_bkg')
    warning('READCDF:noBkg','READCDF: No background flag found');
    FLAG_bkg = logical(zeros(1,n_datapoints));
    COL_bkg = 0;
  end
  % PUT THE CDF DATA INTO A USEABLE STRUCTURE
  for jj=1:n_variables
    if (jj == COL_bkg)
      continue;
    end
%    ID = header.Variables{jj};
    ID = strrep(header.Variables{jj},FieldSeparationString,'.');
    temp = [data_in{:,jj}];
    switch header.Variables{jj,4}
     % EPOCH REQUIRES A SPECIAL CASE
     case 'epoch'
      if ~ISOCTAVE
        % converts CDF date (ms since 01-Jan-0000)
        % to a matlab date (days since 01-Jan-0000)
        temp = todatenum(temp);
      end
      dt_on = temp(~FLAG_bkg);
      dt_off = temp(FLAG_bkg);
      if isempty(dt_on)
        t0 = 0;
      else
        t0 = dt_on(1);
      end
      dataout.dt.data = (dt_on - t0)*TimeConversionFactor;
      dataout.dt.bkg  = (dt_off- t0)*TimeConversionFactor;
      dataout.dt.unit = dt_unit;
      if ISOCTAVE
        oct_t0 = localtime(t0/1000 - 1970*365.2422*86400 - 13.8*3600 + 58);
        oct_datestr = strftime('%d-%b-%Y %T', oct_t0);
        dataout.dt.info = ['\\Delta t from ' oct_datestr];
        disp('Ignore the warning about \\D');
      else
        dataout.dt.info = ['\Delta t from ' datestr(t0)];
      end
      % next line forces dt priority to be same as that of time
      dataout.dt.priority = header.VariableAttributes.priority{jj,2};
      % ALL OTHER DATA TYPES HANDLED HERE
     otherwise
    end
    eval(['dataout.' ID '.data = temp(~FLAG_bkg);']);
    eval(['dataout.' ID '.bkg  = temp(FLAG_bkg);']);
    clear temp;
  end
  % PUT VARIABLE ATTRIBUTES INTO NAMED SUBSTRUCTURES
  var_attr_name = fieldnames(header.VariableAttributes);
  for jj=1:length(var_attr_name)
    %temp is a cell array
    temp = header.VariableAttributes.(var_attr_name{jj});
    for kk=1:length({temp{:,1}})
      if (kk ~= COL_bkg)
	tempVarName = strrep(temp{kk,1},FieldSeparationString,'.');
	eval(['dataout.' tempVarName '.(var_attr_name{jj}) =' ...
		    ' temp{kk,2};' ]);
      end
    end
    clear temp;
  end
  % PUT GLOBAL ATTRIBUTES INTO NAMED SUBSTRUCTURES
  glo_attr_name = fieldnames(header.GlobalAttributes);
  if ISOCTAVE
    for jj=1:length(glo_attr_name)
      glo_attr_struct{jj} = strrep(glo_attr_name{jj},FieldSeparationString,'.');
    end
  else
    glo_attr_struct = strrep(glo_attr_name,FieldSeparationString,'.');
  end
  for jj=1:length(glo_attr_name)
    if ~ischar(header.GlobalAttributes.(glo_attr_name{jj}){1})
      eval(['dataout.globals.' glo_attr_struct{jj} ...
	    '= [header.GlobalAttributes.(glo_attr_name{jj}){:}];']);
    else
      eval(['dataout.globals.' glo_attr_struct{jj} ...
	    '= header.GlobalAttributes.(glo_attr_name{jj});']);
    end
  end
  