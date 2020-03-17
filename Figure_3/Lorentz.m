function y=Lorentz(X,x)
% This function is called by leasqr
% x is a vector which contains the coefficients of the
% equation. X and Y are the option data sets that were
% passed to leasqr.
% A-peak height  B-centre  C-width  D-background
B=x(1);
C=x(2);
A=x(3);
D=x(4);
E = ((X-B).^2)/(C.^2);
y = A./(1+E) + D;
