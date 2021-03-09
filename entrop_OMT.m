%  Created by Sagar K. Tamang (January 2021)
%
%  This program approximates the optimal transportation plan using entropic
%  regularization of the cost function
%  
%  The algorithm used in this code is referenced from the following:
%  Tamang et. al. "Ensemble Riemannian Data Assimilation over the 
%  Wasserstein Space"

function [U] = entrop_OMT(x,y,p,q,gamma,niter)
    x = x'; y = y';
    N(1) = size(x,2);
    N(2) = size(y,2);
    x2 = sum(x.^2,1); y2 = sum(y.^2,1);
    C = repmat(y2,N(1),1)+repmat(x2.',1,N(2))-2*x'*y;
    
    K = exp(-C/gamma);
    
    b = ones(N(2),1); 
    for i=1:niter
        a = p ./ (K*b);
        b = q ./ (K'*a);
    end
    U = diag(a)*K*diag(b);
end
