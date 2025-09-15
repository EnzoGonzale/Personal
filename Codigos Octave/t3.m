function t3
  EA=21*10^(4);
  L=5*10^4;
  pA=8*10^-9;
  ML=1.02;
  dx=10000;
  M=[pA,0,0,0,0;0,pA,0,0,0;0,0,pA,0,0;0,0,0,pA,0;0,0,0,0,ML];
  K=(EA/(dx^2))*[2,-1,0,0,0;-1,2,-1,0,0;0,-1,2,-1,dx/2;0,0,-1,2,-2*dx;0,0,0,-1,(3*dx)/2];
  Ai=inv(inv(M)*K);
  ep=1*10^(-6);
  
  [aval,yn,error,it]=potencia(Ai,ep);
  disp("El menor auto valor es ");disp(aval);
  disp("El auto vector es ");disp(yn);
  disp("El error es ");disp(error);
  disp("El numero de it es ");disp(it);
  Wmenor=(aval)^0.5;
  T=2*pi/Wmenor;
  disp("El W menor es "); disp(Wmenor);
  disp("El periodo T es "); disp(T);
  disp("\n");
  [r,errorr,it]=raices1(ep);
  disp("El valor de lamda es ");disp(r);
  disp("El error es ");disp(errorr);
  disp("El numero de it es ");disp(it);
  Wanalitico=(r/L)*((EA/pA)^0.5);
  disp("El W analitico es "); disp(Wanalitico);
  
endfunction
#Potencia
function [aval,yn,error,it]=potencia(Ai,e)
  y0=[1;1;1;1;1];
  a0=[2;3;1;1;1];#alfa inicial
  it=0;
  itmax=1000;
  fin=0;
  
  while fin==0
    it=it+1;
    xn=y0/max(abs(y0));
    yn=Ai*xn;
    a=yn./xn;#alfa nuevo
    error=max(abs(a-a0)./abs(a));
    if(error<e)
    fin=1;
  endif
  if(it==itmax)
  fin=2;
endif
av2(:,it)=a0;
y0=yn;
a0=a;
aval=1/max(a0);
endwhile
endfunction
#Raizes
function [r,errorr,it]=raices1(error1)
  a=0;
  b=1;
  it=0;
  numsol=0;
  tol=0.2;
  while numsol==0
  it=it+1;
  r=0.5*(a+b);
  errorr=abs(fc(r));
  if (errorr<=error1)
    numsol=1;
  endif
  if it==1000
     numsol=2;
  endif
  if (fc(a)*fc(r)>0)
     a=r;
  else
     b=r;
   endif
 endwhile
endfunction
function f=fc(l)
  f=(1/2547.7)-l*tan(l);
endfunction