function t5
  EA=21*10^4;
  L=5*10^4;
  pA=8*10^(-9);
  px=0;
  ML=1.02*10^(-2);
  dt=0.0002;
  dx=L/5;
  x=0:dx:L;
  my=DifC(L,EA,ML,px,pA);
  [vt,vy1]=runge(L,EA,ML,px,pA);
  #size(vy1)para ver la cantidad de filas y columnas
  MA=zeros(6,5001);
  MA(2:6,:)=vy1;
  A=my-vy1;
  
  for i=1:length(MA(1,:))
  MADX(:,i)=derivadaprimera(MA(:,i),dx);
endfor
for i=1:6
  MADT(i,:)=derivadaprimera(MA(i,:),dt);
endfor


MIT=pA.*(MADT.^2);
for i=1:length(MA)
  IDT(i)=TC(MIT(:,i),dx);
endfor

MIX=EA.*(MADX.^2);
for i=1:length(MA)
  IDX(i)=TC(MIX(:,i),dx);
endfor
X=0.5*IDX;
C=0.5*IDT+0.5*ML*(MADT(6,:).^2);

SX=sum(X);
SC=sum(C);


#Punto 1
plot(vt,my(5,:),"b");
figure;
plot(vt,vy1(5,:),"r");
figure;

#Punto 2
plot(vt,MADT(5,:));
figure;

#Punto 3
plot(vt,C);
figure;

#Punto 4
plot(vt,X);
figure;

#Punto 5
plot(vt,(C+X));
figure;

#Punto 6
plot(vt,A(5,:),"g");
  
disp("Suma de Energia Cinetica: ");disp(SC);
disp("Suma de Energia de Deformacion: ");disp(SX);
endfunction
#Diferencia Central
function my=DifC(L,EA,ML,px,pA)
  tA=0;
  dt=0.0002;
  tf=1;
  dx=L/5;
  M=zeros(5,5);
  for i=1:4
    M(i,i)=pA;
  endfor
  M(5,5)=ML;
  K=(EA/(dx^2))*[2,-1,0,0,0;-1,2,-1,0,0;0,-1,2,-1,0;0,0,-1,2,-1;0,0,dx/2,-2*dx,(3*dx)/2];
  Ua=zeros(5,1);
  dUa=ones(5,1);
  
  
  Uv=Ua-dt*dUa+0.5*(dt^2)*inv(M)*(-K*Ua);#Calculo de y viejo inicial (serie de Taylor)
  it=0;
  
  for i=tA:dt:tf
    it=it+1;
    Un=(dt^2)*inv(M)*(-K*Ua)-Uv+2*Ua;#Ecuacion de Recurrencia (se calcula)
    vx(it)=i;
    my(:,it)=Ua;
    Uv=Ua;
    Ua=Un;
  endfor
endfunction
#Runge Kutta
function [vt,vy1]=runge(L,EA,ML,px,pA)
  t0=0;
  dt=0.0002;
  tf=1;
  dx=L/5;
  y1=zeros(5,1);
  y2=ones(5,1);
  M=zeros(5,5);
  for i=1:4
    M(i,i)=pA;
  endfor
  M(5,5)=ML;
  K=(EA/(dx^2))*[2,-1,0,0,0;-1,2,-1,0,0;0,-1,2,-1,0;0,0,-1,2,-1;0,0,dx/2,-2*dx,(3*dx)/2];
  
  w=0.5;
  it=0;
  
  for i=t0:dt:tf
    it=it+1;
    dy1=y2;
    dy2=inv(M)*(-K*y1);
    k1=dt*dy1;
    k2=dt*dy2;
    
    tg=i+(dt/(2*w));
    y1g=y1+(k1/(2*w));
    y2g=y2+(k2/(2*w));
    dy1g=y2g;
    dy2g=inv(M)*(-K*y1g);
    k1g=dt*dy1g;
    k2g=dt*dy2g;
    
    vt(it)=i;
    vy1(:,it)=y1;
    
    y1=y1+(1-w)*k1+w*k1g;
    y2=y2+(1-w)*k2+w*k2g;
  endfor
endfunction
#Derivada Primera en funcion de t
function VD=derivadaprimera(vy,paso)
  m=length(vy);
  VD(1)=(-3/(2*paso))*vy(1)+(2/(paso))*vy(2)+(-1/(2*paso))*vy(3);
  for i=2:m-1
    VD(i)=(vy(i+1)-vy(i-1))/(2*paso);
  endfor
  VD(m)=(3/(2*paso))*vy(m)+(-2/(paso))*vy(m-1)+(1/(2*paso))*vy(m-2);
endfunction
#Metodo Trapecio
function I=TC(vy,paso)
  m=length(vy);
  w=ones(m,1);
  w(1)=0.5;
  w(m)=0.5;
  w=paso*w;
  I=0;
  for i=1:m
    I=I+w(i)*vy(i);
  endfor
endfunction