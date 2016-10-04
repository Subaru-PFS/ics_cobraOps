function coeff = DerivY_Zern(j, num)

% 

% derivative of Zernikes is expressed as combination of zernikes
% what I am trying to do here is:
%   Given a zernike standard polynomial index j
%   find the corresponding indices of n and m and corresponding angular
%   dependence cosmA or sinmA
%   The derivative of this zernike will be written as sum of individual 
%   zernikes. We need to figure out non-zero zernike terms in this decomposition.
%   we do this by finding out their corresponding np, mp, jp and the
%   coefficients

%   Input:
%       j: positive interger index of zernike polynomials
%       num: number of elements in a vector the coefficients are stored in

% given j, find associated n, m and q. For every jp<j, find the associated 
% np, mp and qp. Determine if the coeeficient is non-zero. Record the jps 
% that have non-zero coefficeints

if(nargin==1)
    num = 50;
end

[n m q] = FindNMQ(j);

coeff = zeros(num,1);
for jp = 1:1:j
    [np mp qp] = FindNMQ(jp);
    if(mod(n-np,2))
        % cal coefficient
        if(abs(m-mp)==1)
            % cal
            if(m*mp==0)
                if ((mod(jp,2)& (m==0)) | (mod(j,2)& (mp==0)))
                coeff(jp) = sqrt(2*(n+1)*(np+1));
                else 
                end
            else
                if (mod(j-jp, 2))
                coeff(jp) = sqrt((n+1)*(np+1));
                    if ((((mp-m)==1)&mod(j,2))|(((mp-m)==-1)& ~mod(j,2)))
                        coeff(jp) = -coeff(jp);
                    else
                    end
                else
                end
            end
        else
        end
    else
    end
end

output = 1:num;
output = output';
output = [output coeff coeff.^2];




