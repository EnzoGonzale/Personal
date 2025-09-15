#Variables
L=0.265; D=0.6; m=0.642; I=0.0343; g=9.81; Lbarra=0.795;
dt=0.01;
t=0:dt:42;
P=zeros(1,length(t));
x0=10*(pi/180); #Angulo Inicial 10 grados

meq=I;
keq=(m*g*D^2)/(4*L);

w=sqrt(keq/meq);
#Fin Variables

#Funciones (Vectores Horizontales)
function x=NewMarckBetaExplicito(dt,t,m,c,k,P,x0,v0,a0)
  x=zeros(1,length(t)); x(1)=x0;
  v=zeros(1,length(t)); v(1)=v0;
  a=zeros(1,length(t)); a(1)=a0;

  kd=k+(3*c)/dt+(6*m)/dt^2;

  for i=1:length(t)-1
    Pd=P(i+1)+m*((6*x(i))/dt^2+(6*v(i))/dt+2*a(i))+c*((3*x(i))/dt+2*v(i)+(dt*a(i))/2);
    x(i+1)=Pd/kd;
    v(i+1)=-2*v(i)-(dt*a(i))/2+(3*(x(i+1)-x(i)))/dt;
    a(i+1)=(1/m)*(P(i+1)-c*v(i+1)-k*x(i+1));
  endfor
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
#Fin Funciones

#Codigo Principal
#Sin Amortiguamiento
x=NewMarckBetaExplicito(dt,t,meq,0,keq,P,x0,0,0); # Tita en funcion del tiempo
y=(Lbarra/2).*sin(x)*100;
#plot(t,y);
#title("Desplazamiento NO Amortiguado");


#Con Amortiguamiento
z=0.00155;
c=2*z*m*w;
ceq=((D^2)/2)*c;

x_amort=NewMarckBetaExplicito(dt,t,meq,ceq,keq,P,x0,0,0);
y_amort=(Lbarra/2).*sin(x_amort)*100;

#figure();
plot(t,y_amort);
title("Desplazamiento Amortiguado");

[y_max, t_max] = BuscarMaximo(t,40,41,y_amort);
y_max
t_max
