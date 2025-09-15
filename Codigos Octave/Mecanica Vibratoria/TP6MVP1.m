function TP6MVP1
m=0.2;
c=0.4;
k=8;
dt=0.12;
t=0:dt:2;
F=zeros(1,length(t));
F(2)=1; F(3)=4; F(4)=9; F(5)=9; F(6)=6;

z=c/(2*(m*k)^0.5);
wn=(k/m)^0.5;
wd=wn*(1-z^2)^0.5;

A=zeros(1,length(t)); A(1)=0;
B=zeros(1,length(t)); B(1)=0;

yc=zeros(1,length(t));
ys=zeros(1,length(t));
x=zeros(1,length(t));

for i=1:length(t)-1
  yc(i)=F(i)*cos(wd*t(i));
  ys(i)=F(i)*sin(wd*t(i));
endfor

for i=1:length(t)-1
  A(i+1)=A(i)*exp(-z*wn*dt)+(dt/(2*m*wd))*(yc(i+1)+yc(i)*exp(-z*wn*dt));
  B(i+1)=B(i)*exp(-z*wn*dt)+(dt/(2*m*wd))*(ys(i+1)+ys(i)*exp(-z*wn*dt));
endfor

for i=1:length(t)-1
  x(i)=A(i)*sin(wd*t(i))-B(i)*cos(wd*t(i));
endfor

Felastica=-k*x

plot(t,x);
endfunction
