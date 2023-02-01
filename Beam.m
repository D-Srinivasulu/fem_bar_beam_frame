E1=3.8e11;Ie1=2.3*10^-8;q0=600;q2=0;F0=10000;M0=4640;L1=0.46;L2=0.38;E2=2.5e11;Ie2=2.4*10^-8;

nele=2;
nnode=nele+1;


ngauss=3;
xivec=[-0.774597 0 0.774597];
wvec=[5/9 8/9 5/9];

coord=[1 0.0;
       2 0.46;
       3 .84];
connect=[1 1 2;
         2 2 3];

P_Load=[2 10000];
P_Moment=[3 4640];

BC_data=[1 1 0;
         1 2 0;
         3 1 0];



K=zeros(2*nnode,2*nnode);
F=zeros(2*nnode,1);

for el=1:nele
    nd1=connect(el,2);
    nd2=connect(el,3);
    vec=[2*nd1-1 2*nd1 2*nd2-1 2*nd2];
    x=[coord(nd1,2) coord(nd2,2)];
    kele=zeros(4,4);
    fele=zeros(4,1);
    for gp=1:ngauss
        
        xi=xivec(gp);w=wvec(gp);
        
        if el==1
            kele(1:4,1:4)=kele(1:4,1:4)+beam_ele_stiff(xi,E1,Ie1,x)*w;
            fele(1:4)=fele(1:4)+beam_ele_load_uniform(xi,q0,L1,x)*w;
        else el==2
            kele(1:4,1:4)=kele(1:4,1:4)+beam_ele_stiff(xi,E2,Ie2,x)*w;
            fele(1:4)=fele(1:4)+beam_ele_load_uniform(xi,0,L2,x)*w;
        
        end
    end
    el
    kele
    fele
    for ii=1:4
        for jj=1:4
            K(vec(ii),vec(jj))=K(vec(ii),vec(jj))+kele(ii,jj);
        end
        F(vec(ii))=F(vec(ii))+fele(ii);
    end
end

for ii=1:size(P_Load,1)
    nd=P_Load(ii,1);
    F0=P_Load(ii,2);
    F(2*nd-1)=F(2*nd-1)+F0;
end

for ii=1:size(P_Moment,1)
    nd=P_Moment(ii,1);
    M0=P_Moment(ii,2);
    F(2*nd)=F(2*nd)+M0;
end

K_glob=K
F_glob=F
supp_dof=[];

for ii=1:size(BC_data,1)
    nd=BC_data(ii,1);
    dof=BC_data(ii,2);
    val=BC_data(ii,3);
    GDOF=2*(nd-1) + dof;
    supp_dof=[supp_dof,GDOF];
    if (val~=0)
        for jj=1:2*nnode
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
K
F
ureduce=inv(K)*F
un=[0;0;ureduce;0]
Freac=K_glob*un

xnume=[];
unume=[];
xi=[-1:0.2:1]';
for el=1:nele
    nd1=connect(el,2);
    nd2=connect(el,3);
    x_n=[coord(nd1,2) coord(nd2,2)];
    vec=[2*nd1-1 2*nd1 2*nd2-1 2*nd2];
    u_n=un(vec);
    
    le=x_n(2)-x_n(1);
    
    Nx=[(1-xi)/2 (1+xi)/2];
    N1=(2-3*xi+xi.^3)/4;
    N2=(1-xi-xi.^2+xi.^3)/4;
    N3=(2+3*xi-xi.^3)/4;
    N4=(-1-xi+xi.^2+xi.^3)/4;
    Nu=[N1 le*N2/2 N3 N4*le/2];
    xnume=[xnume;Nx*x_n'];
    unume=[unume;Nu*u_n];
end   
xnume
unume
