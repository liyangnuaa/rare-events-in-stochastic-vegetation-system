function y=fun2(~,x)

rho=1;
K=10;
beta=3;
x0=1;
alpha=1;
lambda=0.12;
R=1.55;

y=zeros(size(x));
y(1,:)=-rho*x(1,:).*(x(2,:)-x(1,:)/K)+beta*x(1,:)./(x(1,:)+x0);
y(2,:)=-R+alpha*x(2,:)+lambda*x(1,:).*x(2,:);

