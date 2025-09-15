function globalFluidos
#Lista de Funciones:
#Colocar e y D en milimetros y todo en SI
#[Q,v,it]=caso2(D,U,e,L) La Ecuacion de ésta puede variar. Asignarle la formula y variables en el momento
#[D,it]=caso3(U,e,Q) La Ecuacion de ésta puede variar. Asignarle la formula y variables en el momento
#[cm,it]=flujoMasicoIsotermico(P1,P2,D,L,R,T,u,e,f_inicial)
#[P2,it]=presionIsotermica2(P1,P2it,cm,D,L,R,T,u,f)    Al inicio hacer P2it = P1 - 1000
#[P1,it]=presionIsotermica1(P1it,P2,cm,D,L,R,T,u,f)    Al inicio hacer P1it = P2 + 1000
#f=Moody(Re,e,D)  Colocar e y D en metros


f=Moody(40577.3,0.05,77.9)
endfunction





function [Q,v,it]=caso2(D,U,e,L) #U:viscosidad cinematica
#La Ecuacion de ésta puede variar. Asignarle la formula y variables en el momento
tolerancia=0.01;
it=0;
fin=0;
g=9.81;
f_inicial=0.01;
A=(pi*D^2)/4;
er=e/D;
v=0;
while(fin==0)
it=it+1;
#Ecuacion (con f_inicial)
v=((13*2*g)/(1+f_inicial*(L/D)))^(1/2);
#Fin de Ecuacion
Re=(v*D)/U;
f=Moody(Re,e,D);
error=(abs(f-f_inicial)/f);
if(error<=tolerancia)
Q=A*v;
fin=1;
endif
f_inicial=f;
endwhile
endfunction






function [D,it]=caso3(U,e,Q) #U:viscosidad cinematica
#La Ecuacion de ésta puede variar. Asignarle la formula y variables en el momento
tolerancia=0.0001;
it=0;
fin=0;
g=9.81;
f_inicial=0.01;
while(fin==0)
it=it+1;
#Ecuacion (con f_inicial)
D=((f_inicial*Q^2*16)/(pi^2*g*2*0.04))^(1/5);
#Fin de Ecuacion
#Actualizacion de f
Re=(4*Q)/(U*pi*D);
f=Moody(Re,e,D);
error=(abs(f-f_inicial)/f);
if(error<=tolerancia)
fin=1;
endif
f_inicial=f;
endwhile
#Fin Actualizacion
endfunction






function [cm,it]=flujoMasicoIsotermico(P1,P2,D,L,R,T,u,e)
#u:viscosidad, R:constante del gas, Todo en metro
A=(pi*D^2)/4;
f_inicial=0.01;
fin=0;
it=0;
Re=0;
tolerancia=0.0001;
while(fin==0)
it=it+1;
cm = (((P1^2-P2^2)*A^2)/(R*T*(f_inicial*(L/D)+2*log(P1/P2))))^(1/2);
Re=(cm*D)/(u*A);
f=Moody(Re,e,D);
error = (abs(f-f_inicial)/f);
if(error <= tolerancia)
fin=1;
endif
f_inicial=f;
endwhile
endfunction






function [P2,it]=presionIsotermica2(P1,P2it,cm,D,L,R,T,u,f)
#Al inicio hacer P2it = P1 - 1000
it=0;
A=(pi*D^2)/4;
fin = 0;
e = 0.0001;
while (fin == 0)
it=it+1;
P2 = ((P1^2)-(cm^2*R*T/(A^2))*(f*(L/D)+2*log(P1/P2it)))^(1/2);
Error = (abs(P2-P2it)/P2);
P2it = P2;
if ( Error <= e)
  fin = 1;
endif
endwhile
endfunction






function [P1,it]=presionIsotermica1(P1it,P2,cm,D,L,R,T,u,f)
#Al inicio hacer P1it = P2 + 1000
A=(pi*D^2)/4;
fin = 0;
e = 0.0001
it=0;
while (fin == 0)
it=it+1;
P1 = ((P2^2)+(cm^2*R*T/(A^2))*(f*(L/D)+2*log(P1/P2it)))^(1/2);
Error = (abs(P1-P1it))/P1;
P1it = P1;
if ( Error <= e)
  fin = 1;
endif
endwhile
endfunction





function [D,it]=estrechamiento(P1,P2,p,D1,k,cm)
X=P2/P1;
tolerancia=1/100;
fin=0;
it=0;
D2_inicial=0.01;
while (fin==0)
  it=it+1;
  D=((4*cm/pi)*(1-X^(2/k)*(D2_inicial/D1)^4)^(1/2)/(((2*k*P1*p)/(k-1))*(X^(2/k)-X^((k+1)/k)))^(1/2))^(1/2);
  error = (abs(D-D2_inicial)/D);
  if(error<=tolerancia)
  fin=1;
endif
D2_inicial=D;
endwhile
endfunction






function f=Moody(Re,e,D) # Colocar e y D en metros
  er=e/D;
  f = (64./Re).^8 + 9.5.*(log((er./3.7) + (5.74./Re.^0.9)) - (2500./Re).^6).^-16;
  f = f.^0.125;
endfunction
