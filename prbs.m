function u=prbs(n,p,uinit)

% u=prbs(n,p)
% n: shift register length
% p: number of periods in signal

if nargin<3, uinit=ones(1,n);end
    
u=uinit;

q=[0,1,2,1,2,1,1,0,4,3];

if n > 10, error('Maximum allowable value for n is 10');end
if n < 2, error('Minimum allowable value for n is 2');end
if n==8
   for i=n+1:2^n-1, u(i)=xor(xor(u(i-n),u(i-n+1)),xor(u(i-n+2),u(i-n+7)));end
else
    for i=n+1:2^n-1, u(i)=xor(u(i-n),u(i-n+q(n)));end
end

u1=[];
for i=1:p, u1=[u,u1];end

u=2*(u1-0.5)';