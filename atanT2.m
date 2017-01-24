function x = atanT2(y, T, a)
% x = atanT(y,T,a)
% 
% THRESHOLDING FUNCTION OF USING ARCTANGENT PENALTY FUNCTION:
%   gives the solution of 
%   x = argmin_x f(x) = 0.5*(y-x)^2 + T*phi(x,a);
%   where
%   phi(x,a) =  2./(a*sqrt(3))*(atan((2*a*abs(x)+1)/sqrt(3)) - pi/6)
%
% INPUT
%   y : data (scalar or multidimensional array)
%   T : threshold (scalar or multidimensional array)
%   a : penalty convexity parameter (a>0)
%       if a is too small (less than 1e-10) there is no benifit of using the 
%       non-convex penalty function, and the result is approximatly equal
%       to using soft-thresholding.
%
% OUTPUT
%   x : output of atan thresholding

p = soft(y, T);  
x = p;

if ( a >= 1e-10 && T ~= 0 )    

    n = ( p ~= 0 );
    
    yn = y(n);
    absy = abs(yn);
    
    b = 1 - a*absy;
    
    i = ( b == 0 );

    c1 = b.^3./(27.*a^3) - (b.^2)./(6.*a^3) - (absy - T)./(2*a^2) ;
    c2 = b.^2./(9.*a^2) - b./(3.*a^2);
    
    c3 = (sqrt(c1.^2 - c2.^3) - c1).^(1/3);
      
    z = c3 - b/(3*a) + c2./c3;

    if ~isempty(i)
        z(i) = ( (absy(i) - T )/a^2 ) .^ (1/3);
    end

    x(n) = abs(z) .* sign(yn);

end

end

function x = soft(y, T)
% x = soft(y, T)
%
% SOFT THRESHOLDING
% for real or complex data.
%
% INPUT
%   y : data (scalar or multidimensional array)
%   T : threshold (scalar or multidimensional array)
%
% OUTPUT
%   x : output of soft thresholding
%
% If x and T are both multidimensional, then they must be of the same size.

x = max(1 - T./abs(y), 0) .* y;
end