function [x, Tab_a, Tab_k, Tab_r] = algo_OMP(y,H,thres,T)

[N, M] = size(H);

t = 1;
r = y;

x = zeros(1,M);

% Normalizing the columns of H
Htilde = H./sqrt(ones(N,1)*sum(H.^2));

while (norm(r)>thres && t<T)
    
    % Looking for the most correlated column in H with the residual r
    proj = (Htilde'*r);
    [maxVal, ind] = max(abs(proj));
    ind = ind(1);
    
    % Updating the set of non-zero coefficients
    Tab_k(t) = ind;
    hk = H(:,ind);
    
    % Recomputing all non-zero coefficients
    a = (H(:,Tab_k))\y;
    
    % Updating the residual
    r = y-H(:,Tab_k)*a;
    
    % Updating all non-zero components of x
    x(Tab_k) = a;
   
    % saving the evolution of the residual
    Tab_r(t) = norm(r);
    
    % saving the set of non-zero coefficient
    Tab_a{t} = a;

    t=t+1;  
    
end
    
