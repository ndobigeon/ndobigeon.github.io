function [x, Tab_a, Tab_k, Tab_r] = algo_MP(y,H,thres,T)

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
    
    % saving the index of the column
    Tab_k(t) = ind;
    hk = H(:,ind);
    
    % Computing the associated coefficient
    a = hk'*r/norm(hk)^2;
    
    % Updating the residual
    r = r-a*hk;
    
    % Updating the ind-th component of x
    x(ind) = a;
   
    % saving the evolution of the residual
    Tab_r(t) = norm(r);
    
    % saving the set of non-zero coefficient
    Tab_a(t) = a(1);

    t=t+1;

end
    
