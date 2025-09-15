function t4
  EA=21*10^(4);
  L=5*10^4;
  pA=8*10^-9;
  ML=1.02*10^(-2);
  dx=10000;
  M=[pA,0,0,0,0;0,pA,0,0,0;0,0,pA,0,0;0,0,0,pA,0;0,0,0,0,ML];
  K=(EA/(dx^2))*[2,-1,0,0,0;-1,2,-1,0,0;0,-1,2,-1,dx/2;0,0,-1,2,-2*dx;0,0,0,-1,(3*dx)/2];
  Ai=inv(inv(M)*K);
  ep=1*10^(-6);
  
  [aval,yn]=potencia(Ai,ep);
  
  #Punto 1
  W=(aval)^(0.5);
  V=max(yn/max(abs(yn)));
  ti=0;
  tf=10;
  dt=0.1;
  t=ti:dt:tf;
  #Punto 1
  u=V*(cos(W*(t)));
  #Punto 2
  #b
  p1=polnewton(t);
  #c
  tc=ti:dt:1;
  tc2=1:dt:10;
  p2=pollagrange(tc);
  #d
  vx=[0;0.5;1;10];
  vy=[0;1;0;0];
  vx1=[0;0.5]; vy1=[0;1];
  vx2=[0.5;1]; vy2=[1;0];
  vx3=[1;10]; vy3=[0;0];
  td1=0:dt:0.5;
  td2=0.5:dt:1;
  td3=1:dt:10;
  p31=metododirecto(td1,vx1,vy1);
  p32=metododirecto(td2,vx2,vy2);
  p33=metododirecto(td3,vx3,vy3);
  
  #e
  VD=derivadaprimera(p1,dt);
  
  #f
  I=TC(p2,dt);
  disp("Trapecio Compuesto: "); disp(I);
  #g
  IS=SC(p32,dt);
  disp("Simpson Compuesto: "); disp(IS);
  
  #Graficos
  plot(t,u,"b");
  figure;
  plot(t,p1,"r");
  figure;
  hold on
  plot(tc,p2,"g");
  plot(tc2,0,"g");
  hold off
  figure;
  hold on
  plot(td1,p31,"m");
  plot(td2,p32,"m");
  plot(td3,0,"m");
  hold off
  figure;
  plot(t,VD,"b");
endfunction 
#Potencia
function [aval,yn]=potencia(Ai,e)
  y0=[1;1;1;1;1];
  a0=[2;3;1;1;1];#alfa inicial
  it=0;
  itmax=1000;
  fin=0;
  
  while fin==0
    it=it+1;
    xn=y0/max(abs(y0));
    yn=Ai*xn;
    a=yn./xn;#alfa nuevo
    error=max(abs(a-a0)./abs(a));
    if(error<e)
    fin=1;
  endif
  if(it==itmax)
  fin=2;
endif
av2(:,it)=a0;
y0=yn;
a0=a;
aval=1/max(a0);
endwhile
endfunction
#Metodo Newton
function p1=polnewton(x)
  p1=(2*1.*(x-0))+((-4)*1.*(x-0).*(x-0.5))+((8/19)*1.*(x-0).*(x-0.5).*(x-1));
endfunction
#Metodo de Lagrange
function p2=pollagrange(x)
  p2=1*((x-0).*(x-1))/((0.5-0)*(0.5-1));
endfunction
#Metodo Directo
function p3=metododirecto(x,vx,vy)
  n=length(vx);
  for i=1:n
    mf(:,i)=vx.^(i-1);
  endfor
  coef=inv(mf)*vy;
  p3=0;
  for i=1:n
    p3=p3+coef(i)*x.^(i-1);
  endfor
endfunction
#Derivada Primera
function VD=derivadaprimera(vy,paso)
  m=length(vy);
  VD(1)=(-3/(2*paso))*vy(1)+(2/(paso))*vy(2)+(-1/(2*paso))*vy(3);
  for i=2:m-1
    VD(i)=(vy(i+1)-vy(i-1))/(2*paso);
  endfor
  VD(m)=(3/(2*paso))*vy(m)+(-2/(paso))*vy(m-1)+(1/(2*paso))*vy(m-2);
endfunction
#Metodo Simpson
function I=SC(vy,paso)
  m=length(vy);
  w=zeros(m,1);
  w(1)=1/3;
  w(m)=1/3;
  for i=2:2:m-1
    w(i)=4/3;#impar
  endfor
  for i=3:2:m-2
    w(i)=2/3;#par
  endfor
  w=paso*w;
  I=0;
  for i=1:m
    I=I+w(i)*vy(i);
  endfor
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