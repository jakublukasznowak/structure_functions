
function LS = int_ls_short (x,y,maxlag,substract_mean)

if nargin<4
    substract_mean=false;
end

Lx = length(x);
if nargin<3
    maxlag = Lx-1;
end

if nargin<2
    y = x;
end


if substract_mean
    xp = x - mean(x,'omitnan');
    yp = y - mean(y,'omitnan');
else
    xp = x;
    yp = y;
end


Fxy = mean(xp.*yp,'omitnan');


lim = exp(-1);
lag = 0; lag1 = 0; lag2 = 0; 
rho1 = 1; rho2 = 1;

while rho2>lim && lag<maxlag
    
    lag = lag + 1;
    rho = mean( xp(lag+1:Lx) .* yp(1:Lx-lag), 'omitnan' )/Fxy;
    
    if ~isnan(rho)
        lag1 = lag2;
        rho1 = rho2;
        lag2 = lag;
        rho2 = rho;
    end
        
end


LS = interp1( [rho1,rho2], [lag1,lag2], lim );

end