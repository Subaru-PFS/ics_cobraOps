function output=load_dat(filename)
% loads a prelim data set

% $$$   data format is:
% $$$   1  ID
% $$$   2  R.A.          [deg.]
% $$$   3  Dec.          [deg.]
% $$$   4  Exposure Time [sec.]
% $$$   5  Priority      [1(highest)-15(lowest)]
% $$$   6  Magnitude     [AB mag] (currently not available)
% $$$   7  Redshift
% $$$   8  Object Type

   data_format = '%s %f %f %f %d %f %f %s';

   
   fid = fopen(filename);
   contents = textscan(fid, data_format, 'TreatAsEmpty','N/A','CommentStyle','#');

   fclose(fid);
   
   
   output.ra  = contents{2};
   output.dec = contents{3};
   output.exp = contents{4};
   output.pri = contents{5};
   output.mag = contents{6};
   output.z   = contents{7};

   %% just for checking...
   q = output;
   plot((q.ra-mean(q.ra)).*cos(q.dec*pi/180)+mean(q.ra),q.dec, '.');
   xlabel('RA "x cos(Dec)" [deg]');
   ylabel('Dec [deg]');
   title(filename,'Interpreter','none');
   axis equal;