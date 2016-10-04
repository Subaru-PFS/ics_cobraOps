%This function generates the update Vectors according to all that we talked
%about. 
function  [output] = getUpdateVector(qtemp);

%% Get the error/request over angle values back
qtemp;

%% Do the scaling due to error propagation


%% Filter outliers and bad data out


%% et voila, le motor map correction vectors. 

j1fwUv = [];
j1rvUv = [];

j2fwUv = [];
j2rvUv = [];

output.j1fwUv = j1fwUv;
output.j1rvUv = j1rvUv;
output.j2fwUv = j2fwUv;
output.j2rvUv = j2rvUv; 
