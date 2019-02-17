% This function computes the "optimal" correlation lengths and nugget for
% an emulator with design points in the matrix 'Design_points', response
% values in 'y', and predictors in 'H'. The function 'cross_val' is 
% maximised, after this is multiplied (componentwise) with a Gamma density 
% whose parameters are specified as input.
%
% INPUTS
% - Design_points: nx3 matrix with the (orbital) parameters at which the 
%                  values in y have been observed
% - y: nx1 vector of observed outputs
% - H: nxq matrix of predictors (usually has first column of 1s)
% - cor_fun: one of 'exp2', 'matern32', 'matern52', 'abs_exp'
% - m1, m2, m_nu: position of the modes of the gamma distributons used as
%                 prior in the maximisation (respectively for d1, d2, nu).
% - a: shape parameter of the Gamma prior. Smaller a correspond to flatter
%      densities.
%
% OUTPUTS:
% - d:  final value of optimised correlation lengths (2D)
% - nu: final value of optimised nu (1D)


function [d, nu] = max_cross_val(Design_points, y, H, cor_fun, m1, m2, m_nu, a)

M=1; % number of starting points to try
rng(243);
p = haltonset(3,'Skip',floor(1000*rand),'Leap',16); 
% X contains set of starting points for maximisation
X = net(p,M);
X = X*diag(2*[m1,m2,m_nu]); % rescales the columns within a plausible range

%% Part concerning actual maximisation (function h)

% f: Likelihood, approximated through cross validation.
% g: Prior distribution (gamma, with shape parameter alpha>1 and 
%    mode equal to m).
%    The bigger the shape parameter alpha, the more peaked the density.
% h: Main function to be maximised, product of likelihood and prior

f = @(x) cross_val(exp(x(1:2)), exp(x(3)), Design_points, y, H, cor_fun, 'dens');
g = @(x,alpha,m) gampdf(exp(x), alpha, m/(alpha-1));
h = @(x) -sum(log(f(x))) - log( g(x(1),a,m1) * g(x(2),a,m2) * g(x(3),a,m_nu) );

%% Carry out maximisation  starting from two different points
options = optimset('Display', 'off');%, 'TolX', 1e-6, 'TolFun', 1e-8);
x0=log(X(1,:));
[xf,vf] = fminsearch(h, x0, options);

for k=2:M
    x0=log(X(k,:));
    [x_temp,v_temp] = fminsearch(h, log(x0), options);
    if v_temp < vf
        xf=x_temp; vf=v_temp;
    end
end

%% Return the maximisers
d = exp(xf(1:2));
nu = exp(xf(3));

end

