function output = checkLogFile(targetLogFile)

% Read target log
fid  = fopen(targetLogFile,'r');
f=fread(fid,'*char')';
fclose(fid);

% Search log for #INF instances and replace with large number
nInf = length(regexp(f,'#INF'));
if nInf > 0
    disp(sprintf('Found %d instances of #INF. Changing them to 10000',nInf));
    f = regexprep(f,'#INF','10000');
end

% Write the new log text out to the log file
fid  = fopen(targetLogFile,'w');
fprintf(fid,'%s',f);
fclose(fid);

end