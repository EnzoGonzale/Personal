clear; close;
#Variables
r2=0.64; m=4000; c1=c2=2000; k1=k2=20000; L1=1.4; L2=0.9;
I=m*r2;
dt=0.01;
t=0:dt:12;
F=zeros(1,length(t));

M=zeros(2,2); M(1,1)=m; M(2,2)=I;
C=zeros(2,2); C(1,2)=C(2,1)=c2*L2-c1*L1; C(1,1)=c1+c2; C(2,2)=c1*L1^2+c2*L2^2;
K=zeros(2,2); K(1,2)=K(2,1)=k2*L2-k1*L1; K(1,1)=k1+k2; K(2,2)=k1*L1^2+k2*L2^2;

x0=zeros(1,length(M(1,:))); x0(1)=0.05; #El otro termino(x0(2)) es tita=0
v0=zeros(1,length(M(1,:)));
#Fin de Variables

#Funciones (Vectores Horizontales)
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
#Fin de Funciones

#Codigo Principal
[x_t, Xn, X, w, wd, z]=MultiplesGradosLivertadConC(M, K, C, x0, v0, t, dt, F);

plot(t,x_t);
