clear; close;
#Variables
mb=0.642; m=0.00321; mt=mb+m;
L1=0.265; L_barra=0.795; L2=0.2994; #L1=0.265; L2=0.15;
D1=0.795; D2=L_barra;#D1=0.6;
I=0.0343;
g=9.81;
dt=0.01; t=0:dt:100;
F=zeros(1,length(t));

k1=(mt*g)/(2*L1);
k2=(m*g)/L2;


M=zeros(2); M(1,1)=I; M(2,2)=m*(L2^2);
C=zeros(length(M(1,:)));
K=zeros(length(M(1,:))); K(2,1)=K(1,2)=-k2*D2*L2; K(1,1)=0.5*((k1*D1^2)+(k2*D2^2)); K(2,2)=2*k2*(L2^2);

x0=zeros(1,length(M(1,:))); x0(1)=10*(pi/180);
v0=zeros(1,length(M(1,:)));
#Fin Variables

#Funciones
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

function [x_t, Xn, X, w, wd, z]=MultiplesGradosLivertadConC(M, K, C, x0, v0, t, dt, F)
[Xn, lamdas]=eig(K,M);
w=lamdas.^0.5;

Cn=Xn'*C*Xn;
z=zeros(1,length(M(1,:)));
wd=zeros(length(M(1,:)));
for i=1:length(M(1,:))
  z(i)=(1/(2*w(i,i)))*Cn(i,i);
  wd(i,i)=w(i,i)*sqrt(1-z(i));
endfor


for i=0:length(M(1,:))-1
  X(:,i+1)=Xn(:,i+1)/Xn(1,i+1);
endfor

q0=Xn'*M*x0';
qv0=Xn'*M*v0';

for i=1:length(M(1,:))
  for j=0:length(t)-1
    q_t(i,j+1)=exp(-z(i)*w(i,i)*t(j+1))*(cos(wd(i,i)*t(j+1))+(z(i)/(sqrt(1-z(i)^2)))*sin(wd(i,i)*t(j+1)))*q0(i)+(exp(-z(i)*w(i,i)*t(j+1))/wd(i,i))*sin(wd(i,i)*t(j+1))*qv0(i);
  endfor
  [q_CargaGenerica(i,:), yc(i,:), ys(i,:), A(i,:), B(i,:)] = Duhamel(dt,t,w(i,i),wd(i,i),z(i),1,F);
endfor
q_t=q_t+q_CargaGenerica;

x_t=Xn*q_t;

endfunction

function [x_t, Xn, X, w, wd, q_t, q_CargaGenerica]=MultiplesGradosLivertad(z, M, K, x0, v0, t, dt, F)
[Xn, lamdas]=eig(K,M);
w=sqrt(lamdas);
wd=w.*sqrt(1-(z^2));

for i=0:length(M(1,:))-1
  X(:,i+1)=Xn(:,i+1)/Xn(1,i+1);
endfor

q0=Xn'*M*x0';
qv0=Xn'*M*v0';

for i=1:length(M(1,:))
  for j=0:length(t)-1
    q_t(i,j+1)=exp(-z*w(i,i)*t(j+1))*(cos(wd(i,i)*t(j+1))+(z/(sqrt(1-z^2)))*sin(wd(i,i)*t(j+1)))*q0(i)+(exp(-z*w(i,i)*t(j+1))/wd(i,i))*sin(wd(i,i)*t(j+1))*qv0(i);
  endfor
  [q_CargaGenerica(i,:), yc(i,:), ys(i,:), A(i,:), B(i,:)] = Duhamel(dt,t,w(i,i),wd(i,i),z,1,F);
endfor
q_t=q_t+q_CargaGenerica;

x_t=Xn*q_t;

endfunction
#Fin Funciones

#Codigo Principal
x_t=zeros(length(M(1,:)),length(t));
[x_t, Xn, X, w, wd, z]=MultiplesGradosLivertadConC(M, K, C, x0, v0, t, dt, F);
#[x_t, Xn, X, w, wd, q_t, q_CargaGenerica]=MultiplesGradosLivertad(z(1), M, K, x0, v0, t, dt, F);

z=zeros(1,length(M(1,:))); z(1)=0.00155; z(2)=0.002; #Obtenidos Practicamente con decremento logaritmico
c1=2*z(1)*mb*w(1,1); c2=2*z(2)*m*w(2,2);

Ceq=zeros(length(M(1,:))); Ceq(2,1)=Ceq(1,2)=-c2*D2*L2; Ceq(1,1)=0.5*((c1*D1^2)+(c2*D2^2)); Ceq(2,2)=2*c2*L2^2;
Ceq=Ceq+C;

[x_t, Xn, X, w, wd, z]=MultiplesGradosLivertadConC(M, K, Ceq, x0, v0, t, dt, F);

y_t=zeros(1,length(t));

y_t=(L_barra/2).*sin(x_t)*100;

plot(t,y_t(1,:));
title("TMD de dos pendulos simples en los extremos.");
