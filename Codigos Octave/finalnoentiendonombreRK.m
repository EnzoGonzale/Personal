function repaso
  ti=0; dt=1/10; tf=150; td=50;
  t=ti:dt:tf;
  it=0;

  M=zeros(5);
  for i=1:5
    M(i,i)=1;
  endfor
  K=zeros(5);
  for i=1:5
    K(i,i)=2;
  endfor
  for i=1:4
    K(i,i+1)=-1;
    K(i+1,i)=-1;
  endfor
  b=[5;8;9;8;5];
  b=(50/18)*b;

  #Funcion g
  for i=ti:dt:tf
    it=it+1;

    y(it)=g(i);
  endfor
  figure(1);
  plot(t,y);

  #Runge Kutta
  z0=zeros(5,1);
  B=zeros(5,1);
  w=1;
  it2=0;

  for i=ti:dt:tf
    it2=it2+1;

    dz=B;
    dB=b*g(i)-K*z0;
    k1=dt*dz;
    k2=dt*dB;

    tg=i+(dt/(2*w));
    zg=z0+(k1/(2*w));
    Bg=B+(k2/(2*w));
    dzg=Bg;
    dBg=b*g(tg)-K*zg;
    k1g=dt*dzg;
    k2g=dt*dBg;

    zN=z0+(1-w)*k1+w*k1g;
    BN=B+(1-w)*k2+w*k2g;

    m(:,it2)=zN;
    z0=zN;
    B=BN;
  endfor


 #Graficar z3
 z3=m(3,:);
 figure(2);
 plot(t,z3);


 #z en td
 zd=m(:,501); %n=(td/dt)+1
 disp(zd);

 #Funcion discreta u y integral
 xi=0; xf=60; dx=10;
 x=xi:dx:xf;
 it3=1;
 u=zeros(7,1);
 for i=xi:dx:x(5)
   it3=it3+1;
   u(it3)=m(((i/dx)+1),td);
 endfor
 figure(3);
 plot(x,u);

 #Integral
 I=TC((u.^2),dx);
 disp("");
 disp(I);
endfunction
#Funcion g
function y=g(x)
  if x<=50
  y=0*1+(1/25)*1*(x-0)+(-1/625)*1*(x-0)*(x-25);
  else
  y=0;
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
