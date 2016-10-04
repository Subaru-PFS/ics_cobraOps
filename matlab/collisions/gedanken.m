function [err width] = gedanken(alpha,bins)

trials = 1000;

%bins = 100;

dst_per_bin = 2*pi/100;

%alpha = 0.20;
beta = 0.5;

% error per bin-crossing
dx = alpha * dst_per_bin^beta * randn(trials,bins);

% normalized step count needed to reach target
steps = sum(dst_per_bin./(dst_per_bin + dx),2);

err = (mean(steps) - bins)*dst_per_bin;
width = std(steps) * dst_per_bin;
keyboard;
return

disp(format_data(mean(steps)-bins,std(steps),2));
fprintf(1,'alpha * sqrt(bins) = %f\n',alpha * sqrt(bins));


%% things get unstable when alpha*dst_per_bin^(beta-1) becomes
%% "high" (~0.25) so that there is a nonnegligible probability of
%% going backwards.
a = 0.005:0.001:0.10;
for jj =1:length(a)
    ds(jj) = gedanken(a(jj),100);
end
plot(a,ds)
