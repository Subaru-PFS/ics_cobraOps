% starter file for collision investigations

close all;
clear all;

linkLength = 2.375; % arm length
rf = 1.000; % fiber holder radius
patrolRadius = 2*linkLength;
dc = 8.000; % distance between cobras
 
y= sqrt(48);
targets1 = getTargetsAround(0 + 1i* 0);
targets2 = getTargetsAround(8 + 1i* 0);
targets3 = getTargetsAround(-4 + 1i* y);
targets4 = getTargetsAround(4 + 1i* y);
targets5 = getTargetsAround(12 + 1i* y);
targets6 = getTargetsAround(0 + 1i* 2*y);
targets7 = getTargetsAround(8 + 1i* 2*y);

figure(1)
hold on;
plot(targets1,'r.');
plot(targets2,'b.');
plot(targets3,'k.');
plot(targets4,'m.');
plot(targets5,'g.');
plot(targets6,'c.');
plot(targets7,'r.');

%targets1 = 0 + 0.95 *i * rp

targets4tp = XY2TP(targets4-(4 + 1i* y));
targets1tp = XY2TP(targets1);
% figure(4)
% plot(targets1tp,'b.');
elbows4 = targets4 - linkLength * exp(1i* (targets4tp.tht + targets4tp.phi));


elbows1 = targets1 - linkLength * exp(1i* (targets1tp.tht + targets1tp.phi));
range = 3:3;
figure(2)
hold on;
plot([targets4(range);elbows4(range)],'b');
plot(elbows4(range),'bx')
plot(targets4(range), 'b.');
plot(targets1(range), 'r.');
plot(elbows1(range),'rx')
plot([targets1(range);elbows1(range)],'r');
cmplx(@plotcircle,targets1(range),2,'r'); 
cmplx(@plotcircle,(4 + 1i* y),linkLength,'b--');
cmplx(@plotcircle,(4 + 1i* y),patrolRadius,'b--');
cmplx(@plotcircle,(complex(0,0)),linkLength,'r--');
cmplx(@plotcircle,(complex(0,0)),patrolRadius,'r--');
plot(elbows1(1:10), 'r.');


dist1 = pt2linesegment(targets1,elbows4, targets4);
dist2 = pt2linesegment(targets2,elbows4, targets4);
dist3 = pt2linesegment(targets3,elbows4, targets4);
dist5 = pt2linesegment(targets5,elbows4, targets4);
dist6 = pt2linesegment(targets6,elbows4, targets4);
dist7 = pt2linesegment(targets7,elbows4, targets4);

ref = min(dist1,dist2);
ref= min(ref,dist3);
ref = min(ref, dist5);
ref= min(ref,dist6);

disp('Number of ')
collisions = sum(ref<2);
collisions
figure(3)
hold on
plot(targets1(dist1<2),'r.')
plot(targets2(dist2<2),'b.')
plot(targets3(dist3<2),'k.')
plot(targets4(ref<2),'m.')
plot(targets5(dist5<2),'g.')
plot(targets6(dist6<2),'c.')
plot(targets7(dist7<2),'r.')

 %% Test new set of targets using rules: 
 
targets = getTargetsMatrix();
targetsCenterTP = XY2TP(targets(1,:));
elbowsCenter = targets(1,:) - linkLength * exp(1i* (targetsCenterTP.tht + targetsCenterTP.phi));

dist = pt2linesegment(targets(2:end,:), repmat(elbowsCenter,6,1), repmat(targets(1,:),6,1));
collisions = sum(dist<2);
collisions
figure(4)
hold on
plot(targets(dist(1,:)<2),'r.')
 
