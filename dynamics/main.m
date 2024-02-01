clear;
clc;

global rho K beta x0 alpha lambda R sigma1 sigma2 X Y Ux Uy
rho=1;
K=10;
beta=3;
x0=1;
alpha=1;
lambda=0.12;
% R=1.4278;


%% bistable bifurcation curve
NR=1000;
R0=linspace(1.4278,3.5,NR);
node1=zeros(1,NR);
saddle=[];
for i=1:NR
    R=R0(i);
    
    a3=rho*lambda;
    a2=rho*(lambda*x0+alpha);
    a1=rho*alpha*x0+K*beta*lambda-rho*R*K;
    a0=K*beta*alpha-rho*R*K*x0;
    
    p=(3*a3*a1-a2^2)/(3*a3^2);
    q=(27*a3^2*a0-9*a3*a2*a1+2*a2^3)/(27*a3^3);
    
    ommega=(-1+1j*3^0.5)/2;
    x1_1= (-q/2+((q/2)^2+(p/3)^3)^(1/2))^(1/3) + (-q/2-((q/2)^2+(p/3)^3)^(1/2))^(1/3) - a2/(3*a3);
    x1_2= ommega*(-q/2+((q/2)^2+(p/3)^3)^(1/2))^(1/3) + ommega^2*(-q/2-((q/2)^2+(p/3)^3)^(1/2))^(1/3) - a2/(3*a3);
    x1_3= ommega^2*(-q/2+((q/2)^2+(p/3)^3)^(1/2))^(1/3) + ommega*(-q/2-((q/2)^2+(p/3)^3)^(1/2))^(1/3) - a2/(3*a3);
    
    node1(i)=x1_1;
    if real(x1_3)>=0
        saddle=[saddle real(x1_3)];
    end
end
figure;
plot(R0,node1);
hold on
plot(R0(1:length(saddle)),saddle);
hold on
plot([0,R0(length(saddle))],[0,0]);
hold on
plot([1.4278,1.4278],[0,node1(1)]);
hold on
plot([R0(length(saddle)),R0(length(saddle))],[0,node1(length(saddle))]);
hold off

figure;
plot([R0(length(saddle)),4],[0,0]);


%% vector field
R=1.55;

a3=rho*lambda;
a2=rho*(lambda*x0+alpha);
a1=rho*alpha*x0+K*beta*lambda-rho*R*K;
a0=K*beta*alpha-rho*R*K*x0;

p=(3*a3*a1-a2^2)/(3*a3^2);
q=(27*a3^2*a0-9*a3*a2*a1+2*a2^3)/(27*a3^3);

ommega=(-1+1j*3^0.5)/2;
x1_1= (-q/2+((q/2)^2+(p/3)^3)^(1/2))^(1/3) + (-q/2-((q/2)^2+(p/3)^3)^(1/2))^(1/3) - a2/(3*a3);
x1_2= ommega*(-q/2+((q/2)^2+(p/3)^3)^(1/2))^(1/3) + ommega^2*(-q/2-((q/2)^2+(p/3)^3)^(1/2))^(1/3) - a2/(3*a3);
x1_3= ommega^2*(-q/2+((q/2)^2+(p/3)^3)^(1/2))^(1/3) + ommega*(-q/2-((q/2)^2+(p/3)^3)^(1/2))^(1/3) - a2/(3*a3);

x2_1=R/(lambda*x1_1+alpha);
x2_2=R/(lambda*x1_2+alpha);
x2_3=R/(lambda*x1_3+alpha);

SN1=[0;R/alpha];
SN2=[x1_1;x2_1];
US=[real(x1_3);real(x2_3)];

Jacobi=[rho*US(2)-2*rho/K*US(1)-beta*x0/(US(1)+x0)^2, rho*US(1); -lambda*US(2), -alpha-lambda*US(1)];
[Evec,Eval]=eig(Jacobi);
if Eval(1,1)<0
    ev=Evec(:,1);
else
    ev=Evec(:,2);
end
delta=0.01;
Nmf=1000;
X1=[US(1)-ev(1)*linspace(-delta,delta/2,Nmf); US(2)-ev(2)*linspace(-delta,delta/2,Nmf)];
h=0.001;
Nt=6000;
for i=1:Nt
    X2=rk4_2(0,h,X1);
    X1=X2;
end
SM=X2;

if Eval(1,1)>0
    ev=Evec(:,1);
else
    ev=Evec(:,2);
end
X1=[US(1)+ev(1)*linspace(-delta/2,delta,Nmf); US(2)+ev(2)*linspace(-delta/2,delta,Nmf)];
Nt=30000;
for i=1:Nt
    X2=rk4(0,h,X1);
    X1=X2;
end
UM=X2;

figure;
plot(SN1(1),SN1(2),'*');
axis([1 7 0 2])
hold on
plot(SN2(1),SN2(2),'o');
hold on
plot(US(1),US(2),'+');
hold on
plot(SM(1,:),SM(2,:));
hold on
plot(UM(1,:),UM(2,:));
hold off

figure;
fill([7,SM(1,:),7],[4,SM(2,:),0],'m');


%% %%% test
% b1_1=rho*x1_1*(x2_1-x1_1/K)-beta*x1_1/(x1_1+x0);
% b2_1=R-alpha*x2_1-lambda*x1_1*x2_1;

% %%% compute bifurcation value of R         Rc1=1.4278   Rc2=beta*alpha/(rho*x0)=3
% a3=2*rho*lambda;
% a2=rho*(4*lambda*x0+alpha);
% a1=2*rho*(lambda*x0+alpha)*x0;
% a0=rho*alpha*x0^2+K*beta*lambda*x0-K*beta*alpha;
% 
% p=(3*a3*a1-a2^2)/(3*a3^2);
% q=(27*a3^2*a0-9*a3*a2*a1+2*a2^3)/(27*a3^3);
% 
% ommega=(-1+1j*3^0.5)/2;
% x1_1= (-q/2+((q/2)^2+(p/3)^3)^(1/2))^(1/3) + (-q/2-((q/2)^2+(p/3)^3)^(1/2))^(1/3) - a2/(3*a3);
% x1_2= ommega*(-q/2+((q/2)^2+(p/3)^3)^(1/2))^(1/3) + ommega^2*(-q/2-((q/2)^2+(p/3)^3)^(1/2))^(1/3) - a2/(3*a3);
% x1_3= ommega^2*(-q/2+((q/2)^2+(p/3)^3)^(1/2))^(1/3) + ommega*(-q/2-((q/2)^2+(p/3)^3)^(1/2))^(1/3) - a2/(3*a3);
% 
% x2_1=(3*rho*lambda*x1_1^2+2*rho*(lambda*x0+alpha)*x1_1+rho*alpha*x0+K*beta*lambda)/(rho*K);
% x2_2=(3*rho*lambda*x1_2^2+2*rho*(lambda*x0+alpha)*x1_2+rho*alpha*x0+K*beta*lambda)/(rho*K);
% x2_3=(3*rho*lambda*x1_3^2+2*rho*(lambda*x0+alpha)*x1_3+rho*alpha*x0+K*beta*lambda)/(rho*K);


%% collocation points for deep learning
load('SM.mat');
xmin=1;
xmax=7;
ymin=0;
ymax=2;
NH=10000;
x1=(xmax-xmin)*rand(1,NH)+xmin;
x2=(ymax-ymin)*rand(1,NH)+ymin;
J=[];
for i=1:NH
    [m,I]=min(abs(SM(2,:)-x2(i)));
    if x1(i)<SM(1,I)
        J=[J i];
    end
end
x1(J)=[];
x2(J)=[];

figure;
plot(x1,x2,'*');

path = sprintf('x1.mat');
save(path,'x1');
path = sprintf('x2.mat');
save(path,'x2');

N=200;
xlin=linspace(xmin,xmax,N);
ylin=linspace(ymin,ymax,N);
[xmesh,ymesh]=meshgrid(xlin,ylin);
xtest=reshape(xmesh,1,N^2);
ytest=reshape(ymesh,1,N^2);
path = sprintf('xtest.mat');
save(path,'xtest');
path = sprintf('ytest.mat');
save(path,'ytest');


%% computing detH
sigma1=0.1;
sigma2=1;
xnode=SN2;
A=[rho*xnode(2)-2*rho*xnode(1)/K-beta*x0/(xnode(1)+x0)^2,rho*xnode(1);
    -lambda*xnode(2),-alpha-lambda*xnode(1)];
C=zeros(2,2);
C(1,1)=sigma1^2*xnode(1)^4;
C(2,2)=sigma2^2;
A0=[2*A(1,1) 2*A(1,2) 0;A(2,1) A(1,1)+A(2,2) A(1,2);0 2*A(2,1) 2*A(2,2)];
B0=-[C(1,1);C(1,2);C(2,2)];
s=A0\B0;
Z=inv([s(1) s(2);s(2) s(3)]);
detZ=det(Z);


%% NN results
%% loss function
load('LOSS.mat');
figure;
plot(LOSS,'m-');

%% MPEP
load('MPEP.mat');
MPEP=MPEP';
MPEP=[US MPEP SN2];
figure;
plot(MPEP(1,:),MPEP(2,:),'m-');

%% Quasipotential
load('Stest.mat');
xmesh1=reshape(xtest,N,N);
ymesh1=reshape(ytest,N,N);
S=Stest(:,1)';
S=S+(xtest-SN2(1)).^2+(ytest-SN2(2)).^2;
Smesh1=reshape(S,N,N);
load('SM.mat');
for i=1:N
    for j=1:N
        [m,I]=min(abs(SM(2,:)-ylin(j)));
        if xlin(i)<SM(1,I)
            Smesh1(j,i)=inf;
        end
    end
end
lx=Stest(:,2)';
lxmesh1=reshape(lx,N,N);
ly=Stest(:,3)';
lymesh1=reshape(ly,N,N);

figure;
mesh(xmesh1,ymesh1,Smesh1);
figure;
mesh(xmesh1,ymesh1,lxmesh1);
figure;
mesh(xmesh1,ymesh1,lymesh1);

%% characteristic boundary
epsilon_inv=linspace(2,50,1000);
MET_theory=8.4116.*exp(0.1643*epsilon_inv);
ep_inv=[5        10       15       20        25        30        35      40      45];
MET_MC=[10.0548  32.9962  83.4461  192.8937  413.2105  849.3820  1973.8  3842.8  9023.8];
figure;
plot(epsilon_inv,MET_theory);
figure;
plot(ep_inv,MET_MC,'r*');

epsilon_inv2=linspace(2,45,1000);
MET_theory2=10.6968.*exp(0.1763*epsilon_inv2);
ep_inv2=[5        10       15       20        25        30        35      40];
MET_MC2=[12.6107  42.9522  115.9362  237.1639  558.8768  1339.5  2882.0  6570.9];
figure;
plot(epsilon_inv2,MET_theory2);
figure;
plot(ep_inv2,MET_MC2,'r*');

epsilon_inv3=linspace(2,35,1000);
MET_theory3=7.6311.*exp(0.2454*epsilon_inv3);
ep_inv3=[5        10       15       20        25        30];
MET_MC3=[14.5769  60.7287  218.3536  704.8839  1971.8  6145.7];
figure;
plot(epsilon_inv3,MET_theory3);
figure;
plot(ep_inv3,MET_MC3,'r*');

%% non-characteristic boundary
%% Computing prefactor via NN
% epsilon=1/25;
X=xmesh1;
Y=ymesh1;
dx=xlin(2)-xlin(1);
dy=ylin(2)-ylin(1);
H_bar=detZ;
xb=3;
[m,Ix]=min(abs(xlin-xb));
[m,Iy]=min(Smesh1(:,Ix));
x_star=[xlin(Ix);ylin(Iy)];
[Ux,Uy]=gradient(Smesh1,dx,dy);
miu_star=-(0.5*sigma1^2*x_star(1)^4*Ux(Iy,Ix)+lxmesh1(Iy,Ix));
[Uxy,Uyy]=gradient(Uy,dx,dy);
deth_star=Uyy(Iy,Ix);
delta=0.05;
T=10;
Nstep=floor(T/h);
t2=-h;
x00=x_star;
MPEP2=x00;

for i=1:Nstep
    x1=rk4_3(0,h,x00);
    if norm(x1-SN2)<=delta
        break;
    end
    x00=x1;
    MPEP2=[x00 MPEP2];
    t2=[-(i+1)*h,t2];
end
%MPEP=[xnode MPEP];
figure;
plot(MPEP2(1,:),MPEP2(2,:),'m-');

[L1_1,L1_2]=gradient(lxmesh1,dx,dy);
[L2_1,L2_2]=gradient(lymesh1,dx,dy);
Np2=length(t2);
divL2=[];
normbx2=[];
for i=1:Np2
    x=MPEP2(:,i);
    divL2=[interp2(X,Y,L1_1,x(1),x(2),'linear')+interp2(X,Y,L2_2,x(1),x(2),'linear')+4*sigma1^2*x(1)^3*interp2(X,Y,Ux,x(1),x(2),'linear'),divL2];
end
integ2=(sum(divL2))*h;
prefactor2=1/miu_star*sqrt(2*pi*deth_star)/sqrt(abs(det(H_bar)))*exp(integ2);

epsilon_inv=linspace(30,150,1000);
MET2_theory=prefactor2./sqrt(epsilon_inv).*exp(0.0691*epsilon_inv);
ep_inv=[40, 60, 80, 100, 120, 140];
MET2_MC=[123.9558, 482.5652, 1917.3, 7134.3, 30521, 130760];
figure;
plot(epsilon_inv,MET2_theory);
figure;
plot(ep_inv,MET2_MC,'r*');

