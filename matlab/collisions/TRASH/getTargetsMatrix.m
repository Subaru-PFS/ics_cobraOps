function target = getTargetsMatrix()

numtrg = 1000; % Specify number of targets to generate
reachAround = 0.05; % Specify reach factor. Eg: 0.05 means that 5% of
% reachable radial distance is excluded for target generation so that there
% is some padding on ideal keepout zones.
numPos = 7;
L1 = 2.375;
L2 = L1;
phiMax = 180;
phiMin = 0;
% Inner keepout radius
ikor = sqrt(L1^2+L2^2-2*L1*L2*cos(phiMin*pi/180));
% Max reach
reach = sqrt(L1^2+L2^2-2*L1*L2*cos(phiMax*pi/180));

jj = 0:5;

center = repmat((8 * exp(i *jj * pi / 3 ))',1,numtrg);
center = [repmat(complex(0,0),1,numtrg);center];
gamma = rand(numPos,numtrg)*2*pi;
pad = reachAround*reach;
bad  = reshape((1:numtrg*numPos),numtrg,numPos)'; %1:numtrg; % "bad" radii are outside the range
% Riko + pad : reach - pad
radius = repmat((1:numtrg),numPos,1);
radius(1:numtrg)= 0;

%% Rule 1 Inner keep out zone:
while ~isempty(bad)
    radius(bad) = sqrt(rand(size(bad))) .* (reach - pad);
    bad = find(radius < (ikor + pad)); % replace the inner circle ones.
end

%% Rule 2 No targets closer than 2mm
disp('Number of bad targets from rule #2 (>2mm distance)') 
target = center + radius .* exp(1i*gamma);
numbad = [];
while(true) 
   target(1:1,numbad) = complex(0, 0) + sqrt(rand(size(numbad))) .* (reach - pad) .* exp(1i*rand(size(numbad))*2*pi);
     numbad = [];
    for kk = 1:numtrg
        Y = [real(target(1,kk)), imag(target(1,kk))];
        X = [real(target(2:numPos,kk)), imag(target(2:numPos,kk))];
        D = pdist2(X,Y,'euclidean');
        
        if(any(find(D<2)))
            numbad = [numbad, kk];
            % target(1:1,kk) = complex(0, 0) + sqrt(rand(1)) .* (reach - pad) .* exp(1i*rand(1)*2*pi);
        end
    end
    size(numbad)
    if(isempty(numbad))
        break;
    end
end
end