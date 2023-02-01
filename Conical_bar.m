% Problem Parameters
% ------------------
d1=0;d2=0;
L=0;
rho = 7000;
b = 9.81;

g = 9.81;
w=500;
E=70e9;
nu=0.3;
A=(E/(1-nu^2));


% Gauss Points and weights
% -------------------------
xi1 = -0.774597;   w1 = 5/9;
xi2 = 0.774597;    w2 = 5/9;
xi3 = 0;           w3 = 8/9;

% Element1
% --------
xvec = [0.2;0.3;0.4]; xvec1 = xvec;
K1 = stiff(xi1,xvec,d1,d2,L,E)*w1 + stiff(xi2,xvec,d1,d2,L,E)*w2 + stiff(xi3,xvec,d1,d2,L,E)*w3
F1 = loadvec(xi1,xvec,rho,w,A)*w1 + loadvec(xi2,xvec,rho,w,A)*w2 + loadvec(xi3,xvec,rho,w,A)*w3

K=K1(1:2,1:2);
F=F1(1:2);
u=inv(K)*F