function TP5MVP1
k1=600; k2=1200; k3=2400;
m=2;
x0=[0.3,-0.8,0.3];
dt=0.00001;
t=0:dt:1;
t2=0:dt:3;
M=zeros(3,3); M(1,1)=M(2,2)=M(3,3)=m;

K=zeros(3,3);
K = [600, -600, 0; -600, 1800, -1200; 0, -1200, 3600];


[autovec, autoval]=eig(K,M);
disp("Matriz de frecuencias o armonicos:");
w=autoval.^(0.5)
#Normalizacion respecto al primer componente
for i=1:length(x0)
  X(:,i)=autovec(:,i)./autovec(1,i);
endfor
disp("Matriz de auto-vectores normalizada a la primer componente:");
X

#Normalizacion
for i=1:length(x0)
  Xn(:,i)=X(:,i)./(2^(0.5)*norm(X(:,i)));
endfor
disp("Matriz de auto-vectores normalizada:");
Xn

q0=Xn'*M*x0';

for i=1:3
  for j=0:length(t)-1
  q_t(j+1,i)=q0(i)*cos(w(i,i)*t(j+1));
  endfor
endfor
plot(t,q_t);
title("Desplazamientos Generalizados");

x_t=Xn*q_t';
disp("Los desplazamientos en el instante t = 0.53928s sin amortiguamiento:");
x_t(:,53928)
figure();
plot(t,x_t);
title("Vibraciones no amortiguadas");

#Con amortiguamiento
z=0.1;
wd=w.*(1-z^2)^0.5;
for i=1:3
  for j=0:length(t)-1
    q_tamort(j+1,i)=exp(-z*w(i,i)*t(j+1))*(cos(wd(i,i)*t(j+1))+(z/(1-z^2)^0.5)*sin(wd(i,i)*t(j+1)))*q0(i);
  endfor
endfor

x_tamort=Xn*q_tamort';
disp("Los desplazamientos en el instante t = 0.53928s con amortiguamiento:");
x_tamort(:,53928)
figure();
plot(t,x_tamort);
title("Vibraciones Amortiguadas");


#Con carga P
P0=[5000,0,0];
for i=1:3
  for j=0:length(t2)-1
    P(i,j+1)=P0(i)*cos(1.1*w(1,1)*t2(j+1));
  endfor
endfor

Pm=Xn'*P;
Q_t=zeros(length(t2),3);
for i=1:3
  Q_t(:,i)=IntegralDuhamel(z,w(i,i),wd(i,i),Pm(i,:),t2,dt,1);
endfor

for i=1:3
  for j=0:length(t2)-1
    q_tamort2(j+1,i)=exp(-z*w(i,i)*t2(j+1))*(cos(wd(i,i)*t2(j+1))+(z/(1-z^2)^0.5)*sin(wd(i,i)*t2(j+1)))*q0(i);
  endfor
endfor

x_tamortP=Xn*(q_tamort2+Q_t)';
x_tP=Xn*Q_t';  #size(x_tP) 3x30001
x_max=zeros(3,1);
for i=1:3
  for j=2:length(t2)-1
    if(x_tP(i,j) > x_max(i))
    x_max(i)=x_tP(i,j);
    endif
  endfor
endfor
disp("Amplitudes maximas en regimen permanente");
x_max
figure();
plot(t2,x_tamortP);

endfunction

function [x]=IntegralDuhamel(z,wn,wd,P,t,dt,m)

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
endfunction
