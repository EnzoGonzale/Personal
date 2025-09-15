function repaso
  ti=0;
  tf=150;
  dt=1/10;
  t=ti:dt:tf;
  td=75;
  it=0;

  K=zeros(5,5);
  for i=1:5
    K(i,i)=2;
  endfor
  for i=1:4
    K(i,i+1)=-1;
    K(i+1,i)=-1;
  endfor
  b=(50/19)*[5;8;9;8;5];
  M=zeros(5,5);
  for i=1:5
    M(i,i)=5;
  endfor


#Diferencia central
z0=zeros(5,1);
dz0=zeros(5,1);
zv=z0-dt*dz0+0.5*(dt^2)*(inv(M)*b*g(0)-inv(M)*K*z0);

for i=ti:dt:tf
  it=it+1;
  zn=2*z0-zv+(dt^2)*(inv(M)*b*g(i)-inv(M)*K*z0);

  m(:,it)=z0;
  zv=z0;
  z0=zn;
endfor

x1=0:dt:(75/2);
x2=(75/2):dt:75;
x3=75+dt:dt:100;
y1=g(x1);
y2=g(x2);
y3=g(x3);
figure
hold on;
plot(x1,y1,"b");
plot(x2,y2,"b");
plot(x3,y3,"b");
hold off;

figure;
z3=m(3,:);
plot(t,z3);

#Numero de veces q z3(t) cambia de signo
j=0;
for i=td:length(z3)-1
  if z3(i)*z3(i+1)<0
    j=j+1;
  endif
endfor
disp(j);

#Integral
for i=1:751
vecz3(i)=z3(i);
endfor
%I4=TC(vecz3,dt);
xx=0:0.1:1;
figure(99)
plot(xx)
I3=TC2(xx,0.1);

disp("El valor de la integral:"); disp(I3);

#Derivada
VD3=derivadaprimera(z3,dt);
figure;
plot(VD3,z3);
endfunction
#Funcion g
function y=g(i)
  td=75;
  if i>=td
    y=0;
  elseif i>=(td/2)
    y=((2/75)*i)*(-1)+2;
  else
    y=(2/75)*i;
  endif
endfunction
#Integral
function I=TC2(vec,dx)
  m=length(vec);
  w=ones(m,1);
  w(1)=0.5;
  w(m)=0.5;
  for i=1:m
    I=I+dx*vec(i)*w(i);
  endfor
endfunction
#Derivada
function VD=derivadaprimera(vec,dx)
  m=length(vec);
  VD(1)=((-3/2)*vec(1)+(2*vec(2))+(-1/2)*vec(3))/dx;
  for i=2:m-1
    VD(i)=(vec(i+1)-vec(i-1))/(2*dx);
  endfor
  VD(m)=((3/2)*vec(m)+(-2*vec(m-1))+((1/2)*vec(m-2)))/dx;
endfunction
