% 


function [y] = diracdelta_gauss(x,eps)

y = 1/(eps*sqrt(pi))*exp(-1/eps^2*x.^2);


