function [mu,P]=ekf(f,mu,P,h,z,Q,R)

[mu1,F]=jaccsd(f,mu);  
P=F*P*F'+Q;

[z1,H]=jaccsd(h,mu1);

P12=P*H';

K=P12*inv(H*P12+R);
mu=mu1+K*(z-z1);
P=P-K*P12';

function [z,A]=jaccsd(fun,x)
% JACCSD Jacobian through complex step differentiation
% [z J] = jaccsd(f,x)
% z = f(x)
% J = f'(x)

z=fun(x);
n=numel(x);
m=numel(z);
A=zeros(m,n);
h=n*eps;
for k=1:n
    x1=x;
    x1(k)=x1(k)+h*i;
    A(:,k)=imag(fun(x1))/h;
end