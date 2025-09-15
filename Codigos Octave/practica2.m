function repaso
  x0=0;
  y0=2;
  dx=0.005;
  xf=10;
  w=0.5;
  it=0;
  x=x0:dx:xf;
  for i=x0:dx:xf
    it=it+1;
    
    dy=(-0.5*y0)+2*sin(i);
    k1=dx*dy;
    
    xg=i+(dx/(2*w));
    yg=y0+(k1/(2*w));
    dyg=-0.5*yg+2*sin(xg);
    k2=dx*dyg;
    
    vx(it)=i;
    vy(it)=y0;
    
    y0=y0+(1-w)*k1+w*k2;
  endfor
  plot(vx,vy);
  I=SC(vy,dx);
  disp(vy(1));#Primer punto verificar valores con la grafica
  disp(vy(400));#Calcular haciendo 2/(dx+1)
  disp(vy(1000));#5/(dx+1)
  disp(vy(1800));#10/(dx+1)
  disp("");
  disp(I);
endfunction
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