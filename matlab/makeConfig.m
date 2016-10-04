%%ATTENTION NEEDS UPDATING THIS IS OLD

%% CONFIG FILE GENERATOR
clear all

%% TESTBED CONFIGURATION
% Ref positions of fiducials
% CobraConfig.refpos = [147.881491932322 + 1444.39510450085i,1056.60056512391 + 1449.78848207304i,1976.57629634820 + 1456.41031938212i];
CobraConfig.refpos = [1644.1081542969 + 142.8441162109i, 1257.2365722656 + 517.7436523438i, 94.7736587524 + 1642.8116455078i];

% % Pixel Scale
% CobraConfig.Spx = 54.37; % microns/pixel
% 
%Positioner IDs
pId = [1,2,3,4,5];
% 
% % Center coordinates
% centers = [1720.68226342898 + 947.133554238398i,1429.46960720349 + 942.832813111665i,1130.86273329668 + 941.621259318189i,839.209177458729 + 937.628857200814i,551.325818681879 + 936.924354322245i];
% CobraConfig.centers = centers;
% % Stage 2 offsets (arm 1) (px)
% link1s = [44.0397928290352,44.2041405614929,43.6401416415043,47.5348848762267,44.4310781127331];
% CobraConfig.link1s =link1s;
% % Fiber arm lengths (arm 2) (px)
% link2s = [44.5071589023790,39.4023048809253,48.1724176334825,41.3609906747638,44.4610032221656];
% CobraConfig.link2s = link2s;
% % Stage 1 orientation (deg) CW from +X
% s1HOME = [135,180,90,260,345];
% % Stage 1 ROM (deg)
% s1ROM = [378.364862569303,375.313732442161,378.933884218156,369.976691591435,377.810089807251];
% % Stage 1 outer limit
% s1LIM = s1HOME + s1ROM;
% % s1LIM(find(s1LIM>360)) = s1LIM(find(s1LIM>360)) - 360;
% % Stage 2 Outer limit (deg)
% s2LIM = [145.716983137310,159.898239189629,156.686172589452,162.854590192128,169.887413099165];
% % Stage 2 ROM (deg)
% s2ROM = [148.378947120014,137.110472000324,146.232885865083,147.255216449423,158.144870197116];
% % Stage 2 home
% s2HOME = s2LIM - s2ROM;
% 
% M{1} = round(centers*10000)/10000;
% M{2} = round(link1s*10000)/10000;
% M{3} = round(link2s*10000)/10000;
% M{4} = round(s1ROM*10000)/10000;
% M{5} = round(s2LIM*10000)/10000;
% M{6} = round(s1HOME*10000)/10000;
% % M{7} = round(s2LIM*10000)/10000;
% 
% 
% 
% 
% for ii = pId
%     fldID = strcat('pId',num2str(ii));
%     % S1 centers in CCD frame of ref fiducials. This is center of Cobra frame
%     CobraConfig.(fldID).s1Center = centers(ii);
%     % Link1
%     CobraConfig.(fldID).L1 = link1s(ii);
%     % Link2
%     CobraConfig.(fldID).L2 = link2s(ii);
%     % Orientation
%     CobraConfig.(fldID).orientation = s1HOME(ii);
%     % Outer Patrol radius
%     CobraConfig.(fldID).R_ptrl = sqrt(link1s(ii)^2+link2s(ii)^2-2*link1s(ii)*link2s(ii)*cos(s2LIM(ii)*pi/180));
% end
% 
% 
% 
% for ii=2:length(M)
% format short
% M{ii}
% end
% format long g
% M{1}


%% Save config file
save('EMconfig_06-05-14.mat','CobraConfig')