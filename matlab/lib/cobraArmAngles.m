function [theta phi alpha] = cobraArmAngles(P, S1_link, S2_link, left_arm, debug)

% P is complex number representing position of fiber relative to S1 center
% S1 and S2 links are the link lengths of stage 1 and 2 repsectively

if ~exist('debug','var'), debug = false;end;
if ~exist('left_arm','var'), left_arm = false;end;

% Find angle in IMG space CCW from +X
cntrdAngle = angle(P);

% Force all angles to between 0 and 2pi
if cntrdAngle < 0
    cntrdAngle = cntrdAngle + 2*pi;
end

% disp(strcat('matlab angle is ',num2str(cntrdAngle*180/pi)))

% Convert angle to Testbed space CCW from +X
cntrdAngle = 2*pi-cntrdAngle;

% disp(strcat('testbed angle is ',num2str(cntrdAngle*180/pi)))

% Stage 1 arm CW from P since assumed right arm
S1toP = acos((S2_link^2-abs(P)^2-S1_link^2)/(-2*S1_link*abs(P)));
alpha = S1toP;

% Stage 1 arm CCW from +X_axis
% Matlab theta
switch left_arm
    case false
        thetaM = cntrdAngle - S1toP;
        if thetaM < 0
            thetaM = 2*pi+thetaM;
        end
    case true
        thetaM = cntrdAngle + S1toP;
        if thetaM > 2*pi
            thetaM = thetaM - 2*pi;
        end
end

% disp(strcat('testbed theta angle is ',num2str(thetaM*180/pi)))

% MSIM theta (taken CW from +X_axis)
theta = 2*pi-thetaM;
% disp(strcat('msim theta angle is ',num2str(thetaM*180/pi)))

% Stage 2 CW wrt stage 1 arm
phi = acos( (abs(P)^2 - S1_link^2 - S2_link^2 ) / (-2 * S1_link * S2_link) );
if left_arm
    phi = -phi;
end

if debug
    existingFigure('cobraArmAngles')
    circle(0,0,S1_link+S2_link,'k')
    hold on
    plot(P,'rx')
    plot([0,S1_link*cos(theta)],[0,S1_link*sin(theta)],'r')
    text(0,0,num2str(theta*180/pi))
    plot([S1_link*cos(theta),real(P)],[S1_link*sin(theta),imag(P)],'b')
    text(S1_link*cos(theta),S1_link*sin(theta),num2str(phi*180/pi))
    axis equal
end

return