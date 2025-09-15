function repaso
  xi=0; dx=0.05; xf=10;
  x=xi:dx:xf;
  VD=derivadaprimera(fc(x),dx);
  VD2=derivadasegunda(fc(x),dx);
  vec=VD.*VD2;
  
  I=SC(vec,dx);
  disp(I);
  
  plot(x,VD2);
  
  x2=5:dx:8;
  #VD(100)en x=5 VD(160)en x=8 (5/dx+1) y (8/dx+1)
  j=VD(100);
  for i=100:160
    if j<VD(i)
      j=j;
    else
      j=VD(i);
    endif
  endfor
  disp(j);
  
  disp(VD2(4/dx));
    
endfunction
function y=fc(x)
  y=(10./(x+13369))+((x.^2).*cos(2*x)/10);
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