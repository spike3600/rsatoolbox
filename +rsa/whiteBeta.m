function [u_hat,beta_hat,resMS,Sw_raw,shrink]=whiteBeta(Y,X,indices,varargin)

% Estimates beta coefficiencts beta_hat and residuals from raw time series Y
% Estimates the true activity patterns u_hat by applying noise normalization to beta_hat
% INPUT:
%    - Ys        raw timeseries, T by P concatenated across columns for all
%    runs
%    - Xs        design matrix, T by k, note that this is a GLM in which all
%    runs are modelled together
%    - indices is a structure. indices.row contain the row indices and
%    indices.col contain the column indices for each run in the design
%    matrix
% OPTIONS:
%   'normmode':   'runwise': Does the multivariate noise normalisation by run
%                 'overall': Does the multivariate noise normalisation overall
%   'normmethod': 'multivariate': The is the default using ledoit-wolf reg.
%                 'univariate': Performing univariate noise normalisation (t-values)
%                 'none':    No noise normalisation
% OUTPUT:
%    u_hat    estimated true activity patterns (beta_hat after multivariate noise normalization),
%    resMS    residual mean-square - diagonal of the Var-cov matrix, 1 by P
%    Sw_raw   overall voxel error variance-covariance matrix (PxP), before regularisation 
%    beta_hat estimated raw regression coefficients, K*R by P
%    shrink   applied shrinkage factor 
% Alexander Walther, Joern Diedrichsen
% joern.diedrichsen@googlemail.com
% 2/2015, not-SPM dependent, 2/2019
Opt.normmode = 'overall';  % Either runwise or overall
Opt = rsa.getUserOptions(varargin,Opt);

%%% Discard NaN voxels
test=isnan(sum(Y));
if (any(test))
    warning(sprintf('%d of %d voxels contained NaNs -discarding',sum(test),length(test)));
    Y=Y(:,test==0);
end;


[T,Q] = size(X);
partT = nan(T,1);
partQ = nan(Q,1);
Nrun=numel(indices.row);                                     %%% number of runs
for i=1:Nrun
    partT(indices.row{i},1)=i;
    partQ(indices.col{i},1)=i;
end

beta_hat=pinv(X)*Y;                                       %%% ordinary least squares estimate of beta_hat = inv(X'*X)*X'*Y
res=Y-X*beta_hat;                               %%% residuals: res  = Y - X*beta
erdf = size(X,1)-size(X,2); 
switch (Opt.normmode)
    case 'runwise'
        u_hat   = zeros(size(beta_hat));
        % do run-wise noise normalization
        shrink=zeros(Nrun,1);
        for i=1:Nrun
            idxT    = partT==i;             % Time points for this partition 
            idxQ    = partQ==i;             % Regressors for this partition 
            numFilt = 0;   % Number of filter variables for this run 
            [Sw_hat(:,:,i),shrink(i)]=rsa.stat.covdiag(res(idxT,:),sum(idxT)-sum(idxQ)-numFilt-1);                    %%% regularize Sw_hat through optimal shrinkage
            [V,L]=eig(Sw_hat(:,:,i));       % This is overall faster and numerical more stable than Sw_hat.^-1/2
            l=diag(L);
            sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
            u_hat(idxQ,:)=beta_hat(idxQ,:)*sq;
        end
        shrink=mean(shrink);
        
    case 'overall'
        [Sw_hat,shrink]=rsa.stat.covdiag(res,erdf); %%% regularize Sw_hat through optimal shrinkage
        [V,L]=eig(Sw_hat);
        l=diag(L);
        sq = V*bsxfun(@rdivide,V',sqrt(l)); % Slightly faster than sq = V*diag(1./sqrt(l))*V';
        u_hat=beta_hat*sq;
end
if (nargout>1)
    resMS=sum(res.^2)./erdf;
end 
if (nargout>2)
    Sw_raw=res'*res./erdf;
end 
