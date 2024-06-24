function [y]= PGenCDF(x,p)
y=1/2*(1+sign(x).*gamcdf(abs(x).^p/p,1/p,1));
end