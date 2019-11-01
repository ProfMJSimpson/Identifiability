function [L,M,U] = CredibleInterval(X,P)

% Interpolate and smooth
Xi = linspace(min(X),max(X),2000);
Pi = interp1(X,P,Xi,'spline');

[~,index] = max(Pi);
M = Xi(index);

% Find p
for p = linspace(0,max(P),10000)

    if AreaAbove(Xi,Pi,p) < 0.95
        break;
    end
    
end

% Find L and U
L = min(Xi(Pi > p));

U = max(Xi(Pi > p));

end