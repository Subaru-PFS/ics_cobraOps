function output=defineBenchGeometry(centers, useRealMaps, useRealLinks)
% define bench geometry 
% call defineBenchGeometry([],1,1) to read the config file directly.
%
% thetaDIR is > 0 for positive move out of home (opposite of phi)
% thetaDIR is < 0 for negative move out of home (same as phi)
%
% output is a structure with
%    center (Mx1)
%    L1     (Mx1)
%    L2     (Mx1)
%    phiIn  (Mx1)
%    phiOut (Mx1)
%    tht0   (Mx1)
%    NNmap  (MxM), logical
%    rf
%    distCobras
%    minDist
%    pids 
%    mids 
%    S1Nm 
%    S1Pm, S2Pm, S2Nm, binWidth);

KeepOutAngle = 0.1; % radians (~ 5.5 deg)

alpha = 0.07;
beta = 0.50;
pids = [];
mids = [];

% this is an epsilon to make sure that values near zero in theta are
% intepreted on the positive (or negative) side of the cut,
% respectively. Make it negative for same same direction moveouts (positive
% for positive hardstop).
%
% tht1 implemented below is an opposite-sense hard stop.  when
% that's working, theta direction will be specified on a per cobra
% basis, not on a full instrument basis, as is implemented here.
thteps = 2e-10;  %% vesitigial


if useRealMaps | useRealLinks
    CobraConfig = loadCfgXml;
    cx = [];
    cy = [];
    tht0 = []; % hard stop angle for same sense move-out
    tht1 = []; % hard stop angle for opposite sense move-out
    phiIn = [];
    phiOut = [];
    L1 = [];
    L2 = [];
    S1Nm = []; % slow maps
    S2Pm = [];
    S1Pm = [];
    S2Nm = [];
    F1Nm = []; % fast maps
    F2Pm = [];
    F1Pm = [];
    F2Nm = [];
    for ii = 1:length(CobraConfig.ARM_DATA.ARM_DATA_CONTAINER)
        dc = CobraConfig.ARM_DATA.ARM_DATA_CONTAINER{ii};
        pids = [pids; str2num(dc.DATA_HEADER.Positioner_Id.Text)];
        mids = [mids; str2num(dc.DATA_HEADER.Module_Id.Text)];
        cx = [cx; str2num(dc.KINEMATICS.Global_base_pos_x.Text)];
        cy = [cy; str2num(dc.KINEMATICS.Global_base_pos_y.Text)];
        tht0 = [tht0; (str2num(dc.KINEMATICS.CCW_Global_base_ori_z.Text))*pi/180];
        tht1 = [tht1; (str2num(dc.KINEMATICS.CW_Global_base_ori_z.Text))*pi/180];
        phiIn  = [phiIn ; (str2num(dc.KINEMATICS.Joint2_CCW_limit_angle.Text))*pi/180 - pi];
        phiOut = [phiOut; (str2num(dc.KINEMATICS.Joint2_CW_limit_angle.Text)) *pi/180 - pi];
        if isfield(dc.KINEMATICS, 'Link1_Link_Length')
            L1 = [L1; str2num(dc.KINEMATICS.Link1_Link_Length.Text)];
        else
            L1 = [L1; 25];
        end
        if isfield(dc.KINEMATICS, 'Link2_Link_Length')
            L2 = [L2; str2num(dc.KINEMATICS.Link2_Link_Length.Text)];
        else
            L2 = [L2; 25];
        end
        if isfield(dc.SLOW_CALIBRATION_TABLE, 'Joint1_fwd_stepsizes')
            % regions are expected to be standard bin width of 3.6deg
            % phi out of HS in positive direction
            % tht out of HS in negative direction
            values = str2num(dc.SLOW_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text);
            values = 1./values * 3.6;
            S1Pm = [S1Pm; values(3:end)];
            values = str2num(dc.SLOW_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text);
            values = 1./values * 3.6;
            S2Pm = [S2Pm; values(3:end)];
            values = str2num(dc.SLOW_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text);
            values = 1./values * 3.6;
            S1Nm = [S1Nm; values(3:end)];
            values = str2num(dc.SLOW_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text);
            values = 1./values * 3.6;
            S2Nm = [S2Nm; values(3:end)];
        end
        if isfield(dc.FAST_CALIBRATION_TABLE, 'Joint1_fwd_stepsizes')
            % regions are expected to be standard bin width of 3.6deg
            % phi out of HS in positive direction
            % tht out of HS in negative direction
            values = str2num(dc.FAST_CALIBRATION_TABLE.Joint1_fwd_stepsizes.Text);
            values = 1./values * 3.6;
            F1Pm = [F1Pm; values(3:end)];
            values = str2num(dc.FAST_CALIBRATION_TABLE.Joint2_fwd_stepsizes.Text);
            values = 1./values * 3.6;
            F2Pm = [F2Pm; values(3:end)];
            values = str2num(dc.FAST_CALIBRATION_TABLE.Joint1_rev_stepsizes.Text);
            values = 1./values * 3.6;
            F1Nm = [F1Nm; values(3:end)];
            values = str2num(dc.FAST_CALIBRATION_TABLE.Joint2_rev_stepsizes.Text);
            values = 1./values * 3.6;
            F2Nm = [F2Nm; values(3:end)];
        end
    end
    % calculated map ranges (not physical range of motion).  other
    % option is to read it directly from regions.
    map_range.tht = [0 size(S1Pm,2) * 3.6 * pi / 180];
    map_range.phi = [0 size(S2Pm,2) * 3.6 * pi / 180] - pi;
    
    phiOut = phiOut - KeepOutAngle;
    phiIn  = max(phiIn, -pi); % HACK SOLUTION TO BAD PHI HOMES [2016-07-15]
    phiIn  = phiIn + KeepOutAngle;
    center = cx + 1i*cy;
    %% physical parameters
    rf = 1000/90; % fiber holder radius
    distCobras = 8.000 * 1000/90; % distance between cobras for custom geometries.
    minDist = 2 * rf; % minimum distance between fiber and anything
    %% variables defined in this branch that can propagate to
    %% different centers:
    % L1
    % L2
    % phiIn
    % phiOut
    % S[12][FR]m
    configData = packstruct(L1,L2,phiIn,phiOut,S1Nm,S1Pm,S2Pm,S2Nm,distCobras, pids, mids);
else
    map_range.tht = [0 , 112 * 3.6 * pi / 180];
    map_range.phi = [0 , 112 * 3.6 * pi / 180] - pi;
end

if ~isempty(centers)
    %% Cobra center locations
    center = complex(reshape(centers,[],1));
    csize = size(center);
    ONE = ones(csize);
    
    %% physical parameters
    rf = 1.000; % fiber holder radius
    distCobras = 8.000; % distance between cobras for custom geometries.
    minDist = 2 * rf; % minimum distance between fiber and anything

    thtrange = 385 * pi/180; % range of motion of the theta stage.

    tht0 = 2*pi*rand(csize); % same sense hard stop
    tht1 = mod(tht0 + thtrange, 2*pi); % opp sense hard stop
    
    if useRealMaps
        mapAssignment = ceil(rand(length(center),4) * length(configData.L1));
        S1Pm   = configData.S1Pm(mapAssignment(:,1),:);
        S1Nm   = configData.S1Nm(mapAssignment(:,2),:);
        S2Pm   = configData.S2Pm(mapAssignment(:,3),:);
        S2Nm   = configData.S2Nm(mapAssignment(:,4),:);
    else
        S1Pm = ONE*ones(1,112)*50;
        S2Pm = ONE*ones(1,112)*50;
        S1Nm = ONE*ones(1,112)*50;
        S2Nm = ONE*ones(1,112)*50;
    end
    if useRealLinks
        %take L1,L2, phiIn/Out from CobraConfig.
        pix2mm = distCobras/configData.distCobras;
        LasVegas = ceil(rand(length(center),4) * length(configData.L1));
        L1     = configData.L1(LasVegas(:,1)) * pix2mm;
        L2     = configData.L1(LasVegas(:,2)) * pix2mm;
        phiIn  = configData.phiIn(LasVegas(:,3));
        phiOut = configData.phiOut(LasVegas(:,4));
    else 
        linkLength = 2.375;
        L1 = linkLength * ONE;
        L2 = linkLength * ONE;
        
        DeltaPhi = pi - 2*KeepOutAngle; % throw of the phi arm
        phiOut = -KeepOutAngle * (1 + 0.2 * randn(csize));
        phiIn  = phiOut - DeltaPhi;
    end
end

rMin = abs(L1 + L2 .* exp(1i*phiIn));
rMax = abs(L1 + L2 .* exp(1i*phiOut));
binWidth =  2 * pi/100;  %$ this should be defined higher up in the
                         %code.[PHM 20160901

%% parameters for populating the annular patrol areas uniformly
dA     = 1./( (rMax./rMin).^2 - 1 );  % fractional keepout area
rRange = sqrt(rMax.^2 - rMin.^2);

%% by default, phi home is phi in.  existence of mat files *may*
%% change phi home.
phiHome = phiIn;

%% for the realities of the test bench, having phi too far in is
%% inconvenient.  Use a larger value so that the theta arm does not
%% go crazy.

phiHome = phiHome + 0.5;

%% read in phiHome from elsewhere
matfiles = ls2cell('*.mat');
if length(matfiles) > 0
    disp('Warning -- using mat file positions for home in defineBenchGeometry');
    data = loadTestData('',CobraConfig,1); 
    for id = 1:numel(data)
        if ~isempty(data(id).xst)
            kk = find(pids == id);
            phiHome(kk) = mean(data(id).s2.startAngle(1,:)) * pi / 180 - pi;
        end
    end
end

%% home positions (0 = ss, 1 = os)
home0 = center + L1 .* exp(i*tht0) + L2 .* exp(i*(tht0 + phiHome));
home1 = center + L1 .* exp(i*tht1) + L2 .* exp(i*(tht1 + phiHome));

%% theta overlap
tht_overlap = mod(tht1 - tht0 + pi, 2*pi) - pi;

%% generate the nearest neighbor map/structure
epsilon = 1e-9;  % small number to take care of roundoff errors.
CCdist = abs(bsxfun(@minus, center, center.'));
distCent = 2 * median(L1 + L2) * (8.0 / 9.5);
nnMap = sparse(CCdist < 1.5*distCent & CCdist > epsilon);
[row col] = find(nnMap);
%inUpperTriangle = row < col;
NN.row = row;%(inUpperTriangle);
NN.col = col;%(inUpperTriangle);
NN.xy  = mean([center(NN.row), center(NN.col)],2);

%% calculate some parameters for the field geometry
field.cm = mean(center);
field.R  = max(abs(center - field.cm) + rMax);
field.ang = atan(subarray(cmplx(@polyfit,center,1),1));

output = packstruct(center, L1, L2, phiIn, phiOut, tht0, tht1, rMin, rMax, ...
                    dA, rRange,home0,home1,...
                    nnMap, NN, rf, distCobras, minDist, ...
                    S1Nm, S1Pm, S2Pm, S2Nm, map_range, binWidth, ...
                    alpha, beta, ...
                    pids, mids, thteps, tht_overlap, field);
if exist('mapAssignment','var')
    output.map = mapAssignment;
end

end
