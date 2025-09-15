function TP6MVP3
m=0.1;
c=0.2;
k=5;
dt=0.1;
t=0:dt:8;
P=zeros(1,length(t));
P(2)=5; P(3)=8; P(4)=7; P(5)=5; P(6)=3; P(7)=2; P(8)=1;
tolerancia=0.0001;

x=zeros(1,length(t));
v=zeros(1,length(t));
a=zeros(1,length(t));

for i=1:length(t)-1
fin=0;
a(i)=(1/m)*(P(i)-c*v(i)-k*x(i));
it=0;
while fin==0
it=it+1;
vprev=v(i)+(dt/2)*(a(i)+a(i+1));
xprev=x(i)+v(i)*dt+a(i)*((dt^2)/3)+a(i+1)*((dt^2)/6);

a(i+1)=(1/m)*(P(i+1)-c*vprev-k*xprev);

v(i+1)=v(i)+(dt/2)*(a(i)+a(i+1));
x(i+1)=x(i)+v(i)*dt+a(i)*((dt^2)/3)+a(i+1)*((dt^2)/6);

error=max(abs(a(i+1)-a(i))./a(i+1));
if(error<=tolerancia | it>=100)
fin=1;
endif

endwhile

endfor
M=zeros(length(t),3);
M(:,1)=t; M(:,2)=P; M(:,3)=x;
disp("Con NewmarkBeta implicito:");
M(1:5,:)
Felastica=-k*x;
plotyy(t,x,t,Felastica);

xbeta=NewMarckBetaExplicito(dt,t,m,c,k,P,0,0,0);
figure();
plot(t,xbeta);

B(:,1)=t; B(:,2)=P; B(:,3)=xbeta;
disp("Con NewmarkBeta explicito:");
B(1:5,:)
endfunction

function x=NewMarckBetaExplicito(dt,t,m,c,k,P,x0,v0,a0)
  x=zeros(1,length(t)); x(1)=x0;
  v=zeros(1,length(t)); v(1)=v0;
  a=zeros(1,length(t)); a(1)=a0;

  kd=k+(3*c)/dt+(6*m)/dt^2;

  for i=1:length(t)-1
    Pd=P(i+1)+m*((6*x(i))/dt^2+(6*v(i))/dt+2*a(i))+c*((3*x(i))/dt+2*v(i)+(dt*a(i))/2);
    x(i+1)=Pd/kd;
    v(i+1)=-2*v(i)-(dt*a(i))/2+(3*(x(i+1)-x(i)))/dt;
    a(i+1)=(1/m)*(P(i+1)-c*v(i+1)-k*x(i+1));
  endfor
endfunction
