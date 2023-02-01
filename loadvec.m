function  f = loadvec(xi,xvec,rho,w,A)

% This function calculate integrand of element load vector for an element 
% (xvec) at one gauss point (xi)
% =======================================================================

N = [-xi*(1-xi)/2, 1-xi^2, (xi+1)*xi/2];
re = N*xvec;

B = [(xi-0.5), -2*xi, (xi+0.5)];
J = B*xvec; 


f = (J*rho*w^2*(re)/(A))*N';
