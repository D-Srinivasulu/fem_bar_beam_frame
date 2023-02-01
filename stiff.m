function  Ke = stiff(xi,xvec,d1,d2,L,E)

% This function calculate integrand of element stiffness matrix for an 
% element (xvec) at one gauss point (xi)
% =======================================================================
N = [-xi*(1-xi)/2, 1-xi^2, (xi+1)*xi/2];
re = N*xvec;


B = [(xi-0.5), -2*xi, (xi+0.5)];
J = B*xvec;

Ke = (B'*B/(J)^2-N'*B/(re*J)+N'*N/re^2)*J;
