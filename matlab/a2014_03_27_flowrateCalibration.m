%Preprocess data with AutogenerateMotorMaps.m

for ps = [1:length(crm(:,1))]
    pid = sprintf('pId%d',ps) ; 
 intensities.(pid)(intensities.(pid) < 500) = 0;
sums.(pid) = sum(intensities.(pid),2);
end

for ps = [1:length(crm(:,1))]
    pid = sprintf('pId%d',ps) ; 
meansity = mean(sums.(pid)(1:20));
flowRate.(pid) = meansity*2.5/50
end