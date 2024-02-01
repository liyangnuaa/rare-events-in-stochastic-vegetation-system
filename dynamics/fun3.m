function y=fun3(~,x)

global rho K beta x0 alpha lambda R sigma1 sigma2 X Y Ux Uy

bx=[rho*x(1)*(x(2)-x(1)/K)-beta*x(1)/(x(1)+x0);R-alpha*x(2)-lambda*x(1)*x(2)];
y=zeros(size(x));
y(1)=-(bx(1)+sigma1^2*x(1)^4*interp2(X,Y,Ux,x(1),x(2),'linear'))/norm(bx);
y(2)=-(bx(2)+sigma2^2*interp2(X,Y,Uy,x(1),x(2),'linear'))/norm(bx);

