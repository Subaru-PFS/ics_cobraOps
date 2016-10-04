function output = mapFactor(fBE,mapSize)
% determines the map error factor given fractional bin error and
% the size of the map (fibers x bins) based on one estimate of the
% number of steps to get through a bin, per bin.

%fBE = .1; %fractional bin error. realistic is ~0.25
%mapSize = [10000 100]; % 2048x100 for tht

XX = zeros(1,prod(mapSize)); % distance travelled
MF = zeros(1,prod(mapSize)); % map multiplication factor
ctr = 0;

while numel(XX) ~= sum(sum(XX))
    ctr = ctr + 1;
    moveMe = (XX < 1); % logical of map cells to move

    Xremaining = 1 - XX;
    
    dx = inf(size(XX));
    dx(moveMe) = randn(1,sum(moveMe)) * fBE + 1;
    
    Xend = XX + dx; % # bins moved through

    done = Xend >= 1; % bin is done when total travel > 1
    
    fracMove = Xremaining./dx;
    
    fracMove(~done & moveMe) = 1;
    
    MF = MF + fracMove;
    XX = min(Xend,1);
    
    X(ctr,:) = XX; % X holds the move history

end

MF = reshape(MF,mapSize); % ultimately this multiplies against the
                          % inverse motor map (steps/bin)

output = MF;
