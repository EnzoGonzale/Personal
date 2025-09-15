function FINAL_BALDASSO_DAVID

clc;clear;
t0=0;
tf=150;
dt=0.1; 
vt=t0:dt:tf;
z0=[0;0;0;0;0];
dz0=[0;0;0;0;0];
M=5*eye(5);
K=[2,-1,0,0,0;-1,2,-1,0,0;0,-1,2,-1,0;0,0,-1,2,-1;0,0,0,-1,2];
b=(50/19)*[5;8;9;8;5];
%DERIVADA CENTRAL 
zv=z0-dt*dz0+0.5*dt^2*inv(M)*(b*g(t0)-K*z0);
it=0;
for i=t0:dt:tf
  it=it+1;
  zn=(inv(M)*b*g(i)-K*z0)*dt^2-zv+2*z0;
  vz(:,it)=z0;
  zv=z0;
  z0=zn;
endfor 
%GRAFICO DE Z3 EN FUNCION DEL TIEMPO
z3=vz(3,:);
figure(1)
title('La funcion z3 en funcion del tiempo es:')
plot(vt,z3);
%CALCULO DE V3
v3=derivada(z3,dt);
disp(v3);
%GRAFICO DE Z3 EN FUNCION DE V3
figure(2)
title('z3 en duncion de v3');
plot(v3,z3)
%CALCULO DE LA FUNCION C
c=integral(v3,dt)
disp(c);
%GRAFICO DE C EN FUNCION DE V3
figure(3)
title('C en funcion de v3');
plot(v3,c)
%CANTIDAD DE VECES QUE Z3 CAMBIA DE SIGNO
pos1=(75/dt)+1;
pos2=(150/dt)+1;
zz3=z3(pos1:pos2);
cant=0;
n=length(zz3);
for j=2:n
  z3(j-1)*z3(j)<0;
  cant=cant+1;
endfor 
disp('Las veces que la funcion z3 cambia su signo son');
disp(cant);
endfunction
function y=g(t)
  td=75;
  if t<=td/2
    y=t/td;
  endif
  if t>=td/2 && t<=td
    y=-t/td;
  endif
  if t>=td
    y=1;
  endif
endfunction 
function v3=derivada(vector,paso)
  n=length(vector);
  v3(1)=(-3*vector(1))/(2*paso)+(2*vector(2))/paso-vector(3)/(2*paso);
  v3(n)=(3*vector(n))/(2*paso)-(2*vector(n-1))/paso+vector(n-2)/(2*paso);
  for i=2:n-1
    v3(i)=(vector(i+1)-vector(i-1)/2*paso);
  endfor
endfunction
function sum=integral(v3,paso)
  sum=0;
  n=length(v3);
  coef=ones(n,1);
  coef(1)=1/2;
  coef(n)=1/2;
  for i=1:n
    sum=sum+(coef(i)*v3(i));
  endfor
  sum=paso*sum;
endfunction