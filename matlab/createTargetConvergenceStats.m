%% Get the move size for each iteration.
close all;
j1errs = [mId_1_pId_1_str.J1err];

figure(1)
hist(j1errs(1,:), 50)

% Filter to see the distribution

j12 = j1errs(2,j1errs(2,:)<2);
figure(2)
% hist(j12,-0.2:0.01:0.2);
hist(j12,-1:.01:1);
hold on;
j13 = j1errs(3,j1errs(3,:)<2);
hist(j13);
h = findobj(gca,'Type','patch');
%set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')
%set(h,'EdgeColor',[1 0 0],'facealpha',0.5);
set(h(1),'FaceColor',[1 0 0],'EdgeColor',[1 0 0],'facealpha',0.1);

figure(3)
hist(j13,-0.2:.01:0.2)
j14 = j1errs(4,j1errs(4,:)<2);

figure(4)
hist(j14,-0.02:.001:0.02)
j15 = j1errs(5,j1errs(5,:)<2);

figure(5)
hist(j15,-0.02:.001:0.02)

