function TP6MVP2
m=250;
k=4*10^7;
dt=0.001;
t=0:dt:0.3;
P=zeros(1,length(t));
for i=1:length(t)
  if(t(i)<=0.025)
  P(i)=(430*10^3/0.025)*t(i);
  elseif(t(i)<=0.05)
  P(i)=(-430*10^3/0.025)*(t(i)-0.025)+430*10^3;
elseif(t(i)>0.05)
  P(i)=0;
endif
endfor

z=0.05;
wn=(k/m)^0.5;
wd=wn*(1-z^2)^0.5;

A=zeros(1,length(t)); A(1)=0;
B=zeros(1,length(t)); B(1)=0;

yc=zeros(1,length(t));
ys=zeros(1,length(t));
x=zeros(1,length(t));


for i=1:length(t)
  yc(i)=P(i)*cos(wd*t(i));
  ys(i)=P(i)*sin(wd*t(i));
endfor

for i=1:length(t)-1
  A(i+1)=A(i)*exp(-z*wn*dt)+(dt/(2*m*wd))*(yc(i+1)+yc(i)*exp(-z*wn*dt));
  B(i+1)=B(i)*exp(-z*wn*dt)+(dt/(2*m*wd))*(ys(i+1)+ys(i)*exp(-z*wn*dt));
endfor

for i=1:length(t)
  x(i)=A(i)*sin(wd*t(i))-B(i)*cos(wd*t(i));
endfor
plotyy(t,P,t,x);
#hold on
#plot(t,P);
#plot(t,x);
#hold off
endfunction
