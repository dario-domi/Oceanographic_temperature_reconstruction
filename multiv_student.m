% This function generates n samples from the p-dimensional Student distribution.
%
% INPUTS:
% - mu: px1 mean vector
% - C: pxp covariance matrix (the "shape matrix" Sigma satisfies Sigma=(deg-2)/deg * C)
% - deg: degrees of freedom
% - n: size of the generated random sample
%
% OUTPUT:
% - Y: pxn matrix: in each column, an independent sample of the multivariate Student distribution with the above parameters.


function Y = multiv_student(mu, C, deg, n)

p = length(mu);
mu = reshape(mu, [p,1]);

Sigma = (deg-2)/deg * C; % normalisation needed since covariance matrix and kernel matrix are not the same for t-Stud
A = chol(Sigma, 'lower'); %pxp. A*(A') = Sigma

T = mvtrnd(eye(p), deg, n)'; % T is a pxn matrix, with independent p-dimensional samples in each column
% Equivalent to: T = trnd(deg, [p,n]);

% Compute Y by adding mean and a Student sample of right variance
Y = mu + A*T; 

end

