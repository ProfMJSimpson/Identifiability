function A = AreaAbove(varargin)

if nargin == 3

    X = varargin{1};
    Y = varargin{2};
    p = varargin{3};
    
    dX = X(2) - X(1);
 
    A = sum(Y(Y > p)) * dX;
    
elseif nargin == 4
   
    X = varargin{1};
    Y = varargin{2};
    P = varargin{3};
    p = varargin{4};
    
    dA = (X(2) - X(1)) * (Y(2) - Y(1));
    
    A = sum(P(P > p)) * dA;
    
elseif nargin == 5
   
    X = varargin{1};
    Y = varargin{2};
    Z = varargin{3};
    P = varargin{4};
    p = varargin{5};
    
    dA = (X(2) - X(1)) * (Y(2) - Y(1)) * (Z(2) - Z(1));
    
    A = sum(P(P > p)) * dA;
    
end


end