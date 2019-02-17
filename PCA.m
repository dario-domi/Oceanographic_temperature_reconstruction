% This function computes the principal components of a set of n p-dimens vectors (each provided in the form of a N1xN2 
% matrix, where p = N1*N2) and the coefficients of each observation wrt the PCs. Also the StDev of the PCs 
% (as the std of the projections of observations to each PC) and the average input matrix is returned.
% If varargin is present, then it should be a vector of positive weights to be used to modify PCA. (length(weights)=p)
%
% INPUTS
% - X: N1 x N2 x n matrix consisting of n observations. At each N1xN2 level, a matrix of observations.
% - varargin: if present, a vector of positive weights of length p=N1*N2.
%
% OUTPUS
% - PC: N1 x N2 x (n-1) matrix with j-th PC at level j PC(:,:,j)
% - Mn: N1 x N2 matrix, average of initial observations in X
% - Coeff: nx(n-1) matrix. Row Coeff(i,:) stores the coefficients of observation i with respect to the n-1 PCs.
% - Std: (n-1)x1 vector, with standard d canonically associated to the PCs from eigenvalue decomposition.

function [PC, Mn, Coeff, Std] = PCA(X, varargin)

%% General Variables

N1 = size(X,1); N2 = size(X,2); p=N1*N2;
n = size(X,3);

Mn = mean(X,3); % N1 x N2
Xbar = X - Mn; %  N1 x N2 x n
Xbar = reshape(Xbar, [N1*N2, n])'; % n x p
Xbar(isnan(Xbar))=0;

%% PCA: two cases according to whether weights are present or not

if isempty(varargin)  % then do standard PCA
    [U,S,PC] = svd(Xbar, 'econ'); % U: nxn; S: nxn; V: pxn; Y = U*S*V';
    U(:,end) = []; PC(:,end) = []; S(:,end)=[]; S(end,:)=[]; % get rid of data corresponding to last PC
    % U: nx(n-1); S: (n-1)x(n-1); PC: px(n-1); still Xbar - U*S*PC'; 
    Std = diag(S)/sqrt(n-1);
    Coeff = U*S; clear U S

else     % ie, a numeric vector of weights has been provided
    weights=varargin{:};
    if length(weights)~=p
        error('Error: vector of weights has wrong length.');
    end
    
    %% Solving the generalised eigenvalue problem Xbar'*Xbar*W u = \lambda u
    Y=zeros(n,p);
    for i=1:p
        Y(:,i) = Xbar(:,i)*sqrt(weights(i));  % equivalent to Y = Xbar*sqrt(W), but no memory problems
    end
    [U,S,V] = svd(Y, 'econ');  % U: nxn; S: nxn; V: pxn; Y = U*S*V'; 
    U(:,end)=[]; V(:,end)=[]; S(:,end)=[]; S(end,:)=[];
    % Dimensions and relations:
    % U: nx(n-1); S: (n-1)x(n-1); V: px(n-1); still Y = U*S*V'; 
    
    Std = diag(S)/sqrt(n-1);
    for i=1:p
        PC(i,:) = V(i,:) / sqrt(weights(i)); % equivalent to PC=sqrt(W^(-1))*V; 
    end
    clear V;
    Coeff=U*S; clear U S
end

PC = reshape(PC, [N1, N2, n-1]);

end
