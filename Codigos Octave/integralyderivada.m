function integralyderivada
  xi=-2;
  xf=3;
  dx=(xf-xi)/10;
  x=xi:dx:xf;
  x2=xi:1:xf;
  for i=1:length(x)
  y(i)=x(i)^3+x(i)^2-x(i)-1;
endfor
I=SC(y,dx);
disp(I);

for i=1:length(x2)
  y2(i)=x2(i)^3+x2(i)^2-x2(i)-1;
endfor
I2=TC(y2,1);
disp(I2);

VD=derivadaprimera(y,dx)
disp(VD);
disp("VD2");
VD2=derivadasegunda(y,dx);
disp(VD2);
plot(x,VD2,"r");
hold on
plot(x,VD,"b");
plot(x,y,"g");
endfunction
#Metodo Simpson
function I=SC(vy,paso)
  m=length(vy);
  w=zeros(m,1);
  w(1)=1/3;
  w(m)=1/3;
  for i=2:2:m-1
    w(i)=4/3;#impar
  endfor
  for i=3:2:m-2
    w(i)=2/3;#par
  endfor
  w=paso*w;
  I=0;
  for i=1:m
    I=I+w(i)*vy(i);
  endfor
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
#Derivada Segunda
function VD2=derivadasegunda(vy,paso)
  m=length(vy);
  VD2(1)=(2*vy(1)-5*vy(2)+4*vy(3)-1*vy(4))/(paso^2);
  for i=2:m-1
    VD2(i)=(vy(i-1)-2*vy(i)+vy(i+1))/(paso^2);
  endfor
  VD2(m)=(2*vy(m)-5*vy(m-1)+4*vy(m-2)-1*vy(m-3))/(paso^2);
endfunction