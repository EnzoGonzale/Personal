clc; clear; close;

#ParcialMVDhuamel Ej1Resueltos
m=500000;
E=2*10^11;
A=4*10^(-5);
L=1;
F0=100000;
t0=0.4;
z=0.02;

k=(E*A)/L;
wn=sqrt(k/m);
wd=wn*sqrt(1-z^2);

dt=0.001;
t=0:dt:t0;


#FUNCIONES
function [x, yc, ys, A, B]=Duhamel(dt,t,wn,wd,z,m,F)
  yc=zeros(1,length(t));
  ys=zeros(1,length(t));
  A=zeros(1,length(t));
  B=zeros(1,length(t));
  x=zeros(1,length(t));

  for i=1:length(t)
    yc(i)=F(i)*cos(wd*t(i));
    ys(i)=F(i)*sin(wd*t(i));
  endfor

  for i=1:length(t)-1
    A(i+1)=A(i)*exp(-z*wn*dt)+(dt/(2*m*wd))*(yc(i+1)-yc(i)*exp(-z*wn*dt));
    B(i+1)=B(i)*exp(-z*wn*dt)+(dt/(2*m*wd))*(ys(i+1)-ys(i)*exp(-z*wn*dt));
  endfor

  for i=1:length(t)
    x(i)=A(i)*sin(wd*t(i))-B(i)*cos(wd*t(i));
  endfor
endfunction
#FIN DE FUNCIONES


for i=1:length(t)
  if t(i)<=(t0/2)
  F(i)=((F0-0)/((t0/2)-0))*t(i);
  else
  F(i)=((0-F0)/(t0-(t0/2)))*(t(i)-t0);
endif
if t(i)>t0
  F(i)=0;
endif

endfor
plot(t,F);
[x, yc, ys, A, B]=Duhamel(dt,t,wn,wd,z,m,F);
[xnoamort, yc2, ys2, A2, B2]=Duhamel(dt,t,wn,wn,0,m,F);
figure();
plotyy(t,x,t,xnoamort);
xmax=zeros(1,2);
xmaxamort=zeros(1,2);
for i=1:length(t)
  if t(i)<=0.3
    if x(i)>xmax(1)
      xmax(1)=x(i);
      xmax(2)=t(i);
    endif
  endif
  if xnoamort(i)>=xmaxamort(1)
    xmaxamort(1)=xnoamort(i);
    xmaxamort(2)=t(i);
  endif
endfor
xmax
xmaxamort


