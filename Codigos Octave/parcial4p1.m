#Enzo Gonzalez Lugar (1,1)
function parcial4p1
  x1=0; x2=10; dx=0.1;
  x=x1:dx:x2;
  ea=1*10^(-6);

  #Punto1
 [r,it]=raices1(x1,x2,ea);
 disp('La raiz es: ');disp(r);
 
 #Punto2
 disp('La cantidad de iteraciones es: ');disp(it);
 
 #Punto3
 dxn=(x2-r)/(101-1);
 xn=r:dxn:x2;
 I=TC(fcc(xn),dxn);
 disp("");
 disp('El valor de la integral es: ');disp(I);
 
 #Punto4
 x0=-10; dx2=0.4;
 x2=x0:dx2:x2;
 plot(x2,fcc(x2));
endfunction
#Funcion
function fun=fcc(x)
  fun=(x.^2).*cos(x)-x-(1/11300);
endfunction
#Raices
function [r,it]=raices1(a,b,error1)
  it=0;
  itmax=100;
  numsol=0;
  tol=0.2;
  while numsol==0
  it=it+1;
  r=0.5*(a+b);
  errorr=abs(fcc(r));
  if (errorr<=error1)
    numsol=1;
  endif
  if it==itmax
     numsol=2;
  endif
  if (fcc(a)*fcc(r)>0)
     a=r;
  else
     b=r;
   endif
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