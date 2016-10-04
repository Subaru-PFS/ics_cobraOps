%% Define test directory with folders containing the center metrology images
% Currently only looks in subdirectories (need to add in current directory)

% % %% Inputs: Which positioner and filterstring for data
% % pId = 'pId58';
% % FilterString = 'D_10.*SN13';
% % datatitle = 'GEN2 SN13 S1 MM';
% % saveDir = 'C:/Users/cmorantz/Dropbox/PFSMatlab/Results/GEN2_SN13_Center_Metrology/'; %Use / and end with /
% % 
% % % OPTIONAL: Load existing data files that have already been processed by
% % % the centersMain script. NOTE that the data files must currently be in a
% % % subdirectory of this one since the loadPrevious script doesn't look in
% % % dataDir itself, just the subfolders.
% % dataDir = 'C:\Users\cmorantz\Dropbox\PFSMatlab\Results\GEN2_SN13_Center_Metrology'; %Use \ 


if ~exist('loadFilter','var')
    loadFilter = '.';
end

% Get all subdirectories and files
[D F] = subdir(dataDir);

PWD = dir(dataDir);

for file=1:length(PWD)
    if findstr('.mat',PWD(file).name) & findstr(loadFilter,PWD(file).name)
        try
            disp([dataDir '\' PWD(file).name])
            load([dataDir '\' PWD(file).name])
        catch err
            disp(err)
        end
    end 
end


for Dn=1:length(D)
    for Fn=1:length(F{Dn})
        if findstr('.mat',F{Dn}{Fn})>0
            try
                disp([char(D(Dn)) '\' char(F{Dn}{Fn})])
                load([char(D(Dn)) '\' char(F{Dn}{Fn})])
            catch err
                disp(err)
            end
        end
    end
end
    
