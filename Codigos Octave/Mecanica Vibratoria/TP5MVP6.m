#Variables
E=6.9*10^9; I=5.2*10^-6; L=2; m=3000;
k=(3*E*I)/(L^3);
dt=0.001;
t=0:dt:12;
F=zeros(1,length(t));

M=zeros(3); M(1,1)=M(3,3)=m; M(2,2)=4*m;
K=zeros(3); K(1,2)=K(2,1)=K(2,3)=K(3,2)=-k; K(1,1)=K(3,3)=k; K(2,2)=2*k;
x0=zeros(1,length(M(1,:))); x0(1)=0.2;
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

function [x_t, Xn, X, w, wd]=MultiplesGradosLivertad(z, M, K, x0, v0, t, dt, F)
[Xn, lamdas]=eig(K,M);
w=lamdas.^0.5;
wd=w.*sqrt(1-z^2);

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
#Fin de Funciones

#Codigo Principal
[x_t, Xn, X, w, wd]=MultiplesGradosLivertad(0, M, K, x0, v0, t, dt, F);

plot(t,x_t);
