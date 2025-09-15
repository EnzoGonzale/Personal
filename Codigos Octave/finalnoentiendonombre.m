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

  #Diferencia Central
  zA=zeros(5,1);
  dz=zeros(5,1);
  zV=zeros(5,1);
  it2=0;

  for i=ti:dt:tf
    it2=it2+1;

    zN=2*zA-zV+(dt^2)*(b*g(i)-K*zA);
    m(:,it2)=zN;
    zV=zA;
    zA=zN;
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
  y=1*((x-0)*(x-50))/((25-0)*(25-50));
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
