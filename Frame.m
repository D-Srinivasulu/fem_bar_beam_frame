Ie=6.4*10^-6;q0=3200;q2=0;F0=2400;M0=0;L1=0.8;L2=0.6;A=2.4*10^-3;E1=600e6;E2=500e6;E3=100e6;

nele=3;
nnode=nele+1;


ngauss=3;
xivec=[-0.774597 0 0.774597];
wvec=[5/9 8/9 5/9];

coord=[1 0.0 0.0;
       2 0.0 0.8;
       3 0.6 0.8
       4 0.6 0.0];
connect=[1 1 2;
         2 2 3;
         3 3 4];

P_Load_X=[2 2400];
P_Moment=[3 0];

BC_data=[1 1 0;
         1 2 0;
         1 3 0;
         4 1 0;
         4 2 0;
         4 3 0];



K=zeros(3*nnode,3*nnode);
F=zeros(3*nnode,1);

for el=1:nele
    el
    nd1=connect(el,2);
    nd2=connect(el,3);
    vec=[3*nd1-2 3*nd1-1 3*nd1 3*nd2-2 3*nd2-1 3*nd2];
    x=[coord(nd1,2) coord(nd2,2)];
    y=[coord(nd1,3) coord(nd2,3)];
    kele=zeros(6,6);
    fele=zeros(6,1);
    for gp=1:ngauss
        
        xi=xivec(gp);w=wvec(gp);
        if(el==1)
            kele(1:6,1:6)=kele(1:6,1:6)+frame_ele_stiff(xi,E1,Ie,x,y,A)*w;
        elseif el==2
            kele(1:6,1:6)=kele(1:6,1:6)+frame_ele_stiff(xi,E2,Ie,x,y,A)*w;
            fele(1:6)=fele(1:6)+frame_ele_load_uniform(xi,q0,q2,L1,x,y)*w;        
        else
            kele(1:6,1:6)=kele(1:6,1:6)+frame_ele_stiff(xi,E3,Ie,x,y,A)*w;
        end
    end
    
    kele
    fele
    for ii=1:6
        for jj=1:6
            K(vec(ii),vec(jj))=K(vec(ii),vec(jj))+kele(ii,jj);
        end
        F(vec(ii))=F(vec(ii))+fele(ii);
    end
end

for ii=1:size(P_Load_X,1)
    nd=P_Load_X(ii,1);
    F0=P_Load_X(ii,2);
    F(3*nd-2)=F(3*nd-2)+F0;
end
%for ii=1:size(P_Load_Y,1)
%    nd=P_Load(ii,1);
%    F0=P_Load(ii,2);
%    F(3*nd-1)=F(3*nd-2)+F0;
%    end
for ii=1:size(P_Moment,1)
    nd=P_Moment(ii,1);
    M0=P_Moment(ii,2);
    F(3*nd)=F(3*nd)+M0;
end

K_glob=K
F_glob=F
supp_dof=[];

for ii=1:size(BC_data,1)
    nd=BC_data(ii,1);
    dof=BC_data(ii,2);
    val=BC_data(ii,3);
    GDOF=3*(nd-1) + dof;
    supp_dof=[supp_dof,GDOF];
    if (val~=0)
        for jj=1:3*nnode
            F(jj,1)=F(jj,1)-K(jj,GDOF)*val;
        end
    end
end

supp_dof=sort(supp_dof);
for ii=1:size(supp_dof,2)
    dof=supp_dof(ii);
    if dof==1
        K=K(dof+1:end,dof+1:end);
        F=F(dof+1:end);
    elseif dof==2*nnode
        K=K(1:dof-1,1:dof-1);
        F=F(1:dof-1);
    else
        K=K([1:dof-1,dof+1:end],[1:dof-1,dof+1:end]);
        F=F([1:dof-1,dof+1:end]);
    end
    if(ii~=size(supp_dof,2))
        supp_dof(ii+1:end)=supp_dof(ii+1:end)-1;
    end
end

ureduce=inv(K)*F
un=[0;0;0;ureduce;0;0;0]
Freac=K_glob*un

xnume=[];
unume=[];
xi=[-1:0.2:1]';
for el=1:nele
    nd1=connect(el,2);
    nd2=connect(el,3);
    x_n=[coord(nd1,2) coord(nd2,2)];
    y_n=[coord(nd1,3) coord(nd2,3)];
    vec=[3*nd1-2 3*nd1-1 3*nd1 3*nd2-2 3*nd2-1 3*nd2];
    u_n=un(vec);
    
    lex=x_n(2)-x_n(1);
    ley=y_n(2)-y_n(2);
    le=sqrt(lex^2+ley^2);
    Nx=[(1-xi)/2 (1+xi)/2];
    N1=(2-3*xi+xi.^3)/4;
    N2=(1-xi-xi.^2+xi.^3)/4;
    N3=(2+3*xi-xi.^3)/4;
    N4=(-1-xi+xi.^2+xi.^3)/4;
    Nu=[Nx(1) N1 le*N2/2 Nx(2) N3 N4*le/2];
    xnume=[xnume;Nx*x_n'];s
    unume=[unume;Nu*u_n];
end   
xnume
unume
