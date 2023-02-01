function ke=frame_ele_stiff(xi,E,Ie,x,y,A)
lex=x(2)-x(1);
ley=y(2)-y(1);
le=sqrt(lex^2+ley^2);
B1=3*xi/2;
B2=le*(3*xi-1)/4;
B3=-3*xi/2;
B4=le*(3*xi+1)/4;
B=(4/le^2)*[B1 B2 B3 B4];
ke_beam=(E*Ie*le/2)*(B'*B);

N = [(1-xi)/2 (1+xi)/2];
%xe = N*x;
%tht = (d2-d1)/L;
%dx = d1 + tht*xe;
%A = pi*dx^2/4;

Bx = [-1/2 1/2];
len=[sqrt(x(1)^2+y(1)^2) sqrt(x(2)^2+y(2)^2)];
J = Bx*len';

Ke_Bar = A*E*(Bx'*Bx)/J;
ke=zeros(6,6);
l=lex/le
m=ley/le
Q=[l m 0 0 0 0;m -l 0 0 0 0;0 0 1 0 0 0;0 0 0 l m 0;0 0 0 m -l 0;0 0 0 0 0 1];
ke((1:3:4),(1:3:4))=Ke_Bar;
ke([2:3,5:6],[2:3,5:6])=ke_beam;
ke=Q*ke*Q;


