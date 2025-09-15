function repaso
  ti=0; dt=1/1000; tf=4; td=1;
  t=ti:dt:tf;
  it=0;
  pA=8*10^(-9);
  EA=21*10^4;
  Ml=1.02*10^(-2);
  dx=10000; xi=0; xf=50000;
  x=xi:dx:xf;

  M=zeros(5);
  for i=1:4
    M(i,i)=pA;
  endfor
  M(5,5)=-Ml;

  K1=zeros(5);
  for i=1:4
    K1(i,i)=-2;
    K1(i,i+1)=1;
  endfor
  for i=1:3
    K1(i+1,i)=1;
  endfor
  K1=(EA/(dx^2))*K1;

  K2=zeros(5);
  K2(5,3)=0.5; K2(5,4)=-2; K2(5,5)=3/2;
  K2=(EA/dx)*K2;

  K=K1+K2;

  #Funcion g
  for i=ti:dt:tf
    it=it+1;

    y(it)=g(i);
  endfor


  #Diferencia Central
  uA=zeros(5,1);
  uV=zeros(5,1);
  it2=0;

  for i=ti:dt:tf
    it2=it2+1;

    uN=(dt^2)*inv(M)*K*uA+2*uA-uV+(dt^2)*(EA/(dx^2))*10*g(i);
    m(:,it2)=uN;

    uV=uA;
    uA=uN;
  endfor

  n=length(t);
  #u en funcion de t
  u=zeros(6,1);
  for i=1:5
  u(i+1)=m(i,n);
  endfor
  plot(x,u);

  #Calcular matriz dU
  dU=derivadaprimera(u,dt);
  disp(dU);

  #Integral
  vy=(dU.^2)*pA;
  I=TC(vy,dx);
  CIN=(1/2)*I+(1/2)*Ml*(dU(6)^2);

endfunction
function pol=g(x)
  if x<=1
  pol=1*((x-1)*(x-0))/((0.5-1)*(0.5-0));
else
  pol=0;
endif
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
#Derivada Primera
function VD=derivadaprimera(vec,dx)
  n=length(vec);
  VD=zeros(n,1);
  VD(1)=((-3/(2*dx))*vec(1))+(2/dx)*vec(2)-(1/(2*dx))*vec(3);
  for i=2:(n-1)
    VD(i)=(vec(i+1)-vec(i-1))/(2*dx);
  endfor
  VD(n)=(3/(2*dx))*vec(n)-(2/dx)*vec(n-1)+(1/(2*dx))*vec(n-2);
endfunction
