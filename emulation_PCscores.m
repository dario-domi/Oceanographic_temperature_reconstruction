% Function implementing emulation on a dataset defined by Input parameters
% in 'Discr_param' and corresponding outputs in 'Coeff'. 
%
% INPUTS:
% - Design_points: nx3 matrix: in each row, triple of input parameters of 
%      the form: x_emul=(ecc*cos(prec), ecc*sin(prec), obliq) - in radians.
% - Coeff: nxh matrix, with coefficients to emulate in columns
% - r: number of (first) columns of Coeff to actually emulate
% - New_points: Tx3 matrix: in each row, new parameters at which
%               predictions should be made
% - cor_fun: one of the strings 'exp2', 'matern32', 'matern52', 'abs_exp',
%            to specify which correlation function to use
% - var_cov: a string, either 'var' or 'cov'
% - varargin: if present, a matrix rx3. In the first two columns, the corr
%             lengths for the r PCs. In the last column, the nugget term
%
%
% OUTPUTS:
% - M: Txr matrix with emulated means at the paramenter specified in 
%      New_points, wrt the 'input' data in Design_points and Coeff
%
%   Second output changes according to whether 'var' or 'cov' is provided
%   as var_cov input.
%
%   If 'var', then:
% - VarCov: Txr matrix. As M, but with variances rather than means.
%
%   If 'cov', then:
% - VarCov: TxTxr tensor. Cov(:,:,k) represents the covariance matrix of the
%       k-th emulated random vector (built on Design_points and Coeff(:,k))


function [M, VarCov, Dnu] = emulation_PCscores(Design_points, index_lr, Coeff, r, New_points, cor_fun, var_cov, varargin)

%% CHECK THAT SOME INPUTS ARE CORRECT

if ~strcmp(var_cov, 'var') && ~strcmp(var_cov, 'cov')
    error(['Error written by Dario' newline 'Invalid input ''var_cov'':' ...
        'it must be either the string ''var'' or the string ''cov''.']);
end

%%
n = size(Coeff,1);
q = size(Design_points,2) +1; % q-1 regressors will be used
T = size(New_points,1);

%% CHOOSE WHICH REGRESSORS TO USE, AND VALUES OF CORRELATION LENGTHS (IF NOT ALREADY SPECIFIED IN varargin)

Dnu = zeros(r,3);
H_full = cell(r,1); % H_full{c} will contain all discrete regressors for component c
h_full = cell(r,1); % h_full{c} will contain all continuous regressors for component c

%rng(2);
All_discr_regressors = [Design_points, Design_points(:,1:2).^2];
All_cont_regressors =  [New_points, New_points(:,1:2).^2];
for c=1:r
    y=Coeff(:,c);
    
    %% Build matrices of regressors and new parameters
    H = [ones(n,1), All_discr_regressors(:, index_lr(c,:))];
    H_full{c} = H;  % nxq: Matrix of covariates for original data, for component c
    h = [ones(T,1), All_cont_regressors(:, index_lr(c,:))];
    h_full{c} = h;  % Txq: Matrix of covariates for new parameters (ie where interpolation is needed), for component c
    
    %%  Choose correlation lengths
    if isempty(varargin)  % ie, if no corr_lengths have been specified
        m1=0.01;
        m2=0.002;
        nu_initial=0.5;
        a=4;
        [d, nu]=max_cross_val(Design_points, y, H, cor_fun, m1, m2, nu_initial, a);
        disp(['In cross validation, component number ' num2str(c) ' done.']);
        Dnu(c,1:2) = d;
        Dnu(c,3)=nu;
    else  % ie, a specific d was given in input, in varargin
        Dnu = varargin{:};
    end
end


%% REAL EMULATION

M=zeros(T,r);
if strcmp(var_cov, 'var')
    VarCov=zeros(T,r);
else
    VarCov=zeros(T,T,r);
end

for c=1:r
    
    d = [Dnu(c,1), Dnu(c,1), Dnu(c,2)];
    nu = Dnu(c,3);
    y = Coeff(:,c);
    A = Corr_fun(Design_points, Design_points, d, nu, cor_fun);    % nxn
    H = H_full{c};  h = h_full{c};
    K = H'/A;               % qxn, K = H'*(A^-1)
    B = K*H;                % qxq, B = H'*(A^-1)*H
    b = B\(K*y);            % qx1, b = (B^-1)*Ky
    f = y - H*b;            % nx1
    e = A\f;                % nx1, e = A^-1 (y - Hb)
    s2 = (f'*e)/(n-q-2);    % scalar
    
    t = Corr_fun(New_points, Design_points, d, 0, cor_fun);     % Txn
    M(:,c) = (h*b) + (t*e);
    
    p = h' - K*t';          % qxT, p = h(x) - H'*A^-1*t
    
    if strcmp(var_cov, 'var')
        v1 = Corr_fun(Design_points(1,:), Design_points(1,:), d, nu, cor_fun);
        VarCov(:,c) = s2*(v1 - diag(t*(A\t')) + diag(p'*(B\p)) );
    else
        v1 = Corr_fun(New_points, New_points, d, nu, cor_fun);
        VarCov(:,:,c) = s2*(v1 - t*(A\t') + p'*(B\p) );
    end
    
end

end