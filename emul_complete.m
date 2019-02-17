% Given the output of 'emulation_PCscores.m', this function computes temperature emulated mean and covariance for a set
% of Nloc locations. Values of the PCs at the N locations are provided in PC_Val. The average of the original simulations 
% at these locations is provided in Mn_Val.
% Either a full covariance or just the variance at the times corresponding to the rows of M_Pc is output, according to
% whether 'cov' or 'var' is specified as last input.
% The variance from the unused components is added to the variance arising from the linear combination.
%
% INPUTS:
% M_Pc:    Txn matrix (n>=r). In j-th column, emulated mean of the coefficient of the j-th PC.
% Cov_Pc:  Usually, TxTxr. Cov(:,:,j) represents the covariance matrix of the coefficient of the j-th PC. 
%          (last dimension can be >r)
%          If input var_cov = 'var', then Cov_Pc should be a Txr matrix, with only emulated variances in each column.
% PC_Val:  Nloc x n matrix. In each row, values of the n PCs at the location corresponding to that row.
% Mn_Val:  Nloc x 1, with average value of the original simulation outputs at each of the N locations. This will be added 
%          to the linear combination of the columns of Pc_Val.
% Std_PCA: Row or column vector of length n, with standard deviation of PCs (from eigenvalue decomposition).
% r:       number of PCs to use in linear combinations
% var_cov: a string, either 'var' or 'cov'. Accordingly, only variance or full covariance will be computed as output Cov.
%
%
% OUTPUTS:
% M: T x Nloc, with emulated mean temperature at time i and location j stored in M(i,j).
%
% Second output changes according to whether 'var' or 'cov' is provided as var_cov input.
%
% If 'var', then:
%     VarCov: T x Nloc. As M, but with variances rather than means.
%
% If 'cov', then:
%     VarCov: T x T x Nloc. At level (:,:,j), the posterior covariance
%             matrix at j-th location.


function [M, VarCov] = emul_complete(M_Pc, Cov_Pc, PC_Val, Mn_Val, Std_PCA, r, var_cov)

%% Checking last input and coherence with Variance/Covariance provided

if ~strcmp(var_cov, 'var') && ~strcmp(var_cov, 'cov')
    error(['Error written by Dario' newline 'Invalid input ''var_cov'':' ...
        'it must be either the string ''var'' or the string ''cov''.']);
end

if strcmp(var_cov, 'var') && length(size(Cov_Pc))~=2
    error(['Error written by Dario' newline ...
           'Since last input is string ''var'', only variances should' ...
           ' be provided as second input (T x Nloc matrix).']);
end

if strcmp(var_cov, 'cov') && length(size(Cov_Pc))~=3
    error(['Error written by Dario' newline ...
           'Since last input is string ''cov'', full covariance should' ...
           ' be provided as second input (T x T x Nloc matrix).']);
end


% General variables and reshaping

n = length(Std_PCA);   % total number of PCs
T = size(M_Pc,1);      % length of vector containing times
Nloc = size(PC_Val,1); % total number of locations

Std_PCA = reshape(Std_PCA, [1,n]);
M_Pc = M_Pc(:, 1:r);   % Txr


% Building full covariance (or just variance) matrix of PCs, considering "constant" variance for non-used PCs

if strcmp(var_cov, 'var')
    Full_PC_Var = zeros(T,n);    % Txn
    Full_PC_Var(:, 1:r) = Cov_Pc(:, 1:r);
    if r<n
        Full_PC_Var(:,r+1:n) = ones(T,1)* (Std_PCA(r+1:n).^2);
    end    
else
    Full_PC_Cov = zeros(T,T,n);      % TxTxn
    Full_PC_Cov(:, :, 1:r) = Cov_Pc(:,:,1:r);
    if r<n
        for k = r+1:n
            Full_PC_Cov(:, :, k) = 1.e-5 * (Std_PCA(k)^2)* sparse(eye(T));
        end
    end
    
end


% Computing Emulator mean and variances

M = ones(T,1)*Mn_Val' + (M_Pc*PC_Val(:, 1:r)'); % TxNloc + (Txr)x(rxNloc)
Squared_PC = (PC_Val').^2;     % n x Nloc
if strcmp(var_cov, 'var')
    VarCov = Full_PC_Var*Squared_PC; % T x Nloc, from (Txn)x(nxNloc)
else
    VarCov = zeros(T,T,Nloc);
    for k=1:Nloc
        VarCov(:,:,k) = multiprod(Full_PC_Cov, Squared_PC(:,k), [0 3], [1 0]);
    end
    % Equivalent to the following multiprod: TxTxNloc, from (TxTxn)x(nxNloc)
    % VarCov = squeeze(multiprod(Full_PC_Cov, Squared_PC, [0 3], [1 0]));
    % Multiprod however needs to store the matrix of dim T x T x n x Nloc
end

end

