function repaso
ti=0; dt=10; td=750; tf=6*td;
t=ti:dt:tf;
it=0;
for i=ti:dt:1000
  it=it+1;
  y(it)=g(i);
endfor

K=zeros(6);
for i=1:6
  K(i,i)=-2;
endfor
for i=1:5
  K(i,i+1)=1;
  K(i+1,i)=1;
endfor
K(6,5)=2/3; K(6,6)=-2/3;
K=(1/50)*K;

p=zeros(6,1); p(1)=1;
p=(1/50)*p;

#RungeK
w=1;
z0=zeros(6,1);
it2=0;

for i=ti:dt:tf
  it2=it2+1;

  dz=K*z0+p*g(i);
  k1=dt*dz;

  tg=i+(dt/(2*w));
  zg=z0+(k1/(2*w));
  dzg=K*zg+p*g(i);
  k2=dt*dzg;

  zn=z0+(1-w)*k1+w*k2;

  m(:,it2)=zn;

  z0=zn;
endfor
#Graficar z3
z3=m(3,:);
plot(t,z3);

#Vector z(tk)
tk=(tf/2);
tk=(tk/dt)+1;
zj=m(:,tk);
disp(zj);

#Armar y Graficar funcion discreta
x0=1; dx=0.1; xf=1.7;
x=x0:dx:xf;
u=zeros(8,1);
it3=0;

u(1)=g(tf/2);
for i=2:7
  u(i)=m(i-1,tk);
endfor
u(8)=(4*m(6,tk)-m(5,tk))/3;
du=derivadaprimera(u,dx);
figure(2);
plot(x,u,"r");
figure(3);
plot(x,du,"b");

#Integral
du2=du.^2;
I=TC(du2,dx);
disp(I);

#Aproximacion
Fi=ones(8,3);
it4=0;
for i=x0:dx:xf
  it4=it4+1;

  Fi(it4,2)=1/i;
  Fi(it4,3)=1/(i^2);
endfor
a=inv(Fi'*Fi)*(Fi'*u);
it5=0;
for i=x0:dx:xf
  it5=it5+1;

  aprox(it5)=a(1)*1+a(2)*f1(i)+a(3)*f2(i);
endfor
disp(a);
disp("");
disp(aprox);

figure(4);
plot(x,aprox);


endfunction
#Funcion g
function g=g(x)
    if x<750
    g=(2/5)*x;
  else
    g=300;
  endif
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
function y1=f1(x)
  y1=1/x;
endfunction
function y2=f2(x)
  y2=1/(x^2);
endfunction
