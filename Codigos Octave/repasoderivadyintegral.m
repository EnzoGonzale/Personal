function repasoderivadyintegral
  xi=0; dx=10/50; xf=10;
  ti=0; dt=2/100; tf=2;
  x=xi:dx:xf;
  t=ti:dt:tf;
  m=zeros(51,101);
  for i=1:length(x)
    for j=1:length(t)
      m(i,j)=sin(2*x(i)-5*t(j));
    endfor
  endfor
  #Punto 3
  y3=m(6,51:101);
  I3=TC(y3,dt);
  disp(I3);
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
endfunction