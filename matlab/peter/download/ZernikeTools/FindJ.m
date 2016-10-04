function j = FindJ(n, m)

% for a given pair of n and m, find the possible corresponding Js
if n<m
    error('n (1st argument) must be greater than m (2nd argument)!');
end

if mod(n-m, 2)
    error('n (1st argument) minus m (2nd argument) must be even!');
end

% if n is even, there are n+1 terms. if n is odd, there are n+1 terms as
% well

% total number of terms with n' < n
num1 = n*(n+1)/2;

% total number of terms with n' = n but m' < m
if m==0
    num2 = 0;
    j = num1 + 1;
else
    num2 = m-1;
    j = [num1+num2+1, num1+num2+2];
end


