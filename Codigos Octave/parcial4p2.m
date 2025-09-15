#Enzo Gonzalez Lugar (1,1)
function parcial4p2
  xi=0; dx=2/50; xf=2;
  ti=0; dt=10/4; tf=10;
  x=xi:dx:xf;
  t=ti:dt:tf;
  m=zeros(51,5);
  for i=1:length(x)
    for j=1:length(t)
      m(i,j)=5*sin(7*x(i)+0.5*t(j));
    endfor
  endfor
  #Punto1
  for i=1:5
    VDx(:,i)=derivadaprimera(m(:,i),dx);
  endfor
  size(VDx)
  vec1=t.*VDx(16,:); #16=(0.6/dx)+1
  I1=TC(vec1,dt);
  disp("");
  disp("I del punto 1:");
  disp(I1);
  
  #Punto2
  for i=1:51
    VDt(i,:)=derivadaprimera(m(i,:),dt);
  endfor
  
  vec2=x'.*VDt(:,3); #3=(5/dt)+1
  I2=TC(vec2,dx);
  disp("");
  disp("I del punto 2:");
  disp(I2);
  
  #Punto3
  vec3=m.^2;
  for i=1:5
  I3(:,i)=TC(vec3(:,i),dx);
  endfor
  disp("");
  disp("I del punto 3:");
  disp(I3);
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
function VD=derivadaprimera(vy,paso)
  m=length(vy);
  VD(1)=(-3/(2*paso))*vy(1)+(2/(paso))*vy(2)+(-1/(2*paso))*vy(3);
  for i=2:m-1
    VD(i)=(vy(i+1)-vy(i-1))/(2*paso);
  endfor
  VD(m)=(3/(2*paso))*vy(m)+(-2/(paso))*vy(m-1)+(1/(2*paso))*vy(m-2);
endfunction