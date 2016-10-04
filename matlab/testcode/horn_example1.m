[x y] = meshgrid(1:5, 1:1);
ref   = reshape(x + i*y, [], 1);

n_test = 10000;
sigma_test = 0.1; % gaussian error to introduce

errs = sigma_test * (randn(length(ref), n_test) + ...
                     i*randn(length(ref), n_test));

R = randn(1,n_test) * .001;
S = randn(1,n_test) * .001 + 1;
T = (randn(1,n_test) + i*randn(1,n_test)) * sigma_test * .01;

fids = bsxfun(@plus, ref, errs);

%apply rotation and scale
%fids = bsxfun(@times, fids, S .* exp(i*R));
%apply translation
%fids = bsxfun(@plus, fids, T);


for jj = 1:n_test
    hr(jj) = horn87(fids(:,jj), ref);
    h1(jj) = horn87(fids(:,jj), fids(:,1));
end

R1 = [h1.R];
S1 = [h1.S];
T1 = [h1.T];

fixed = bsxfun(@times, fids, S1 .* exp(i*R1));
fixed = bsxfun(@plus, fixed, T1);

Rr = [hr.R];
Sr = [hr.S];
Tr = [hr.T];

bootref = mean(fids,2);
