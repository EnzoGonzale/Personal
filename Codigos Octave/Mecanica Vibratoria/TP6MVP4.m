function TP6MVP4
m=0.2;
c=0.4;
k=8;
dt=0.12;
t=0:dt:6;
dt2=0.006;
t2=0:dt2:6;
F=zeros(1,length(t));
F(2)=1; F(3)=4; F(4)=9; F(5)=9; F(6)=6;
P=interp1(t,F,t2,'linear');
z=c/(2*(m*k)^0.5);
wn=(k/m)^0.5;
wd=wn*(1-z^2)^0.5;

x=zeros(1,length(t2));
v=zeros(1,length(t2));
a=zeros(1,length(t2));

for i=1:length(t2)-1
  x(i+1)=x(i)+dt2*v(i)+((dt2^2)/(2*m))*(P(i)-c*v(i)-k*x(i));
  v(i+1)=(2/dt2)*(x(i+1)-x(i))-v(i);
endfor
x(1:6)
Felastica=-k*x;
#plot(t2,P);
plotyy(t2,x,t2,Felastica);
endfunction

