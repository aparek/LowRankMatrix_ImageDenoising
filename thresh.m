function [ y ] = thresh( x, lam, a, pen, p )
% y = thresh(x,lam)
% Generalized thresholding in case of BP with log/atan penalty
% Ankit Parekh

y = zeros(size(x));
ind = abs(x)>=lam;

if nargin == 4
    p = 0;
end
if ~lam
    y = x;
else

    switch pen
        case 'firm'
            y = ( 1/(1-a*lam)*max(abs(x)-lam,0) + (1-1/(1-a*lam))*max(abs(x)-1/a,0) ).*sign(x);
        case 'log'
            y(ind) = (abs(x(ind))/2-1/(2*a) + sqrt((abs(x(ind))./2 + 1/(2*a)).^2 - lam/a)).*sign(x(ind));
        case 'atan'
            y = atanT2(x,lam,a);
        case 'l1'
            y = soft(x,lam);
        case 'lp'
            y(ind) = max((abs(x(ind)) - lam.^(2-p).*abs(x(ind)).^(p-1)),0).*sign(x(ind));
        case 'TS1'
            y = TS1thresh(x,lam,a);
        otherwise
            disp('Please select penalty from the following: log, atan, l1')
    end
end


