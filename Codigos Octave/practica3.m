function repaso
  x1=5; x2=20; x0=0;
  x=x0:0.1:x2;
  ea=0.00001;

  xr=sec(x1,x2,ea);
 disp('La raiz es: ');disp(xr);
 xn=x0:0.1:xr;
 I=TC(fcc(xn),0.1);
 disp("");
 disp(I);

 xr2=sec2(3,5,ea);
 disp(fcc(xr2));
 hold on
 plot(x,fcc(x));
 hold off
endfunction
#Funcion
function fun=fcc(x)
  fun=(x.^2).*cos(x)-x-(1/13369);
endfunction
#Derivada
function fun2=dfcc(x)
  fun2=(2*x).*cos(x)-(x^2).*sin(x)-1;
endfunction
function xr=sec(x1,x2,ea)
  it=0;
  itmax=1000;
  numsol=0;
  while numsol==0
  it=it+1;
  m=(fcc(x2)-fcc(x1))/(x2-x1);
  xr=x2-(fcc(x2)/m);
  errorr=abs(fcc(xr));
  if (errorr<=ea)
    numsol=1;
  endif
  if it==itmax
     numsol=2;
  endif
  x1=x2;
  x2=xr;
  e(it)=errorr; #error como vector
  xrr(it)=xr; #raiz como vector
 endwhile
endfunction
function xr2=sec2(x1,x2,ea)
  it=0;
  itmax=1000;
  numsol=0;
  while numsol==0
  it=it+1;
  m=(dfcc(x2)-dfcc(x1))/(x2-x1);
  xr2=x2-(dfcc(x2)/m);
  errorr=abs(dfcc(xr2));
  if (errorr<=ea)
    numsol=1;
  endif
  x1=x2;
  x2=xr2;
 endwhile
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