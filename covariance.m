function [C] = covariance(X,Y)
    N = size(X,2);
    X_bar = mean(X,2);
    Y_bar = mean(Y,2);
    
    
    C = (X-X_bar*ones(1,N))*(Y-Y_bar*ones(1,N))';
    C = C./(N-1);
end
