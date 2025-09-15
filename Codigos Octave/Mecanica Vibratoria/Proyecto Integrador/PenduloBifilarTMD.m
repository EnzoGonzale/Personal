#Variables
m1=0.642; m2=0.353; mt=m1+m2;
L1=0.265; L2=0.15; L_barra_grande=0.795; L_barra_chica=0.39; #L1=0.265; L2=0.15;
D1=0.6; D2=0.13; #D1=0.6; D2=0.21;
I1=0.0343; I2=0.0045;
g=9.81;
dt=0.01; t=0:dt:10;
F=zeros(1,length(t));
#zitta=0,0295 para barra a 20,5 cm


Meq=zeros(2,2); Meq(1,1)=I1; Meq(2,2)=I2; #Meq(1,2)=Meq(2,1)=-2*I2;
Keq=zeros(2,2); Keq(1,2)=Keq(2,1)=-(m2/L2)*D1*D2; Keq(1,1)=((mt/L1)+(m2/L2))*D1^2; Keq(2,2)=(m2/L2)*D2^2;
Keq=(g/4)*Keq;

x0=zeros(1,length(Meq(1,:))); x0(1)=x0(2)=10*(pi/180); #Angulos iniciales de 10 grados
v0=zeros(1,length(Meq(1,:)));
#Fin Variables

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

function [x_max, t_max] = BuscarMaximo(t,ti,tf,x)
  x_max=0;
  t_max=0;
  i_inicial=0; fin_inicial=0;
  i_final=0; fin_final=0;
  for j=1:length(t)-1
    if t(j)>=ti && fin_inicial==0
      i_inicial=j;
      fin_inicial=1;
    endif
    if t(j)>=tf && fin_final==0
      i_final=j;
      fin_final=1;
    endif
  endfor

  for i=i_inicial:i_final
    if x(i)>x_max
      x_max=x(i);
      t_max=t(i);
    endif
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
#Fin Funciones

#Codigo Principal
#Sin Amortiguamiento
[x_t, Xn, X, w, wd, q_t, q_CargaGenerica]=MultiplesGradosLivertad(0, Meq, Keq, x0, v0, t, dt, F);

y_t=zeros(2,length(t));
y_t(1,:)=(L_barra_grande/2).*sin(x_t(1,:))*100;
y_t(2,:)=(L_barra_chica/2).*sin(x_t(2,:))*100;

#plot(t,y_t);
#title("Desplazamiento NO Amortiguado");


#Con Amortiguamiento
z=zeros(1,2);
z(2)=0.0295;
z(1)=0.00155;
c1=2*z(1)*m1*w(1,1);
c2=2*z(2)*m2*w(2,2);
Ceq=zeros(length(Meq(1,:))); Ceq(2,1)=Ceq(1,2)=-c2*(D1*D2)/2; Ceq(1,1)=(c1+c2)*((D1^2)/2); Ceq(2,2)=c2*((D2^2)/2);
z_nuevoo

[x_t_amort, Xn, X, w, wd, z_nuevoo]=MultiplesGradosLivertadConC(Meq, Keq, Ceq, x0, v0, t, dt, F);

y_t_amort=zeros(2,length(t));
y_t_amort(1,:)=(L_barra_grande/2).*sin(x_t_amort(1,:))*100;
y_t_amort(2,:)=(L_barra_chica/2).*sin(x_t_amort(2,:))*100;

#figure();
plot(t,y_t_amort);
title("Desplazamiento Amortiguado con TMD Modificado", 'fontsize', 18);
set(gca, 'fontsize', 16);

