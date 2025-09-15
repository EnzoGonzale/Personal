function actividad3
  clc
  
  ra=8*10^-9;
  M=[ra,0,0,0,0;0,ra,0,0,0;0,0,ra,0,0;0,0,0,ra,0;0,0,0,0,1.02];
  ea=21*10^4;
  l=5*10^4;
  dx=10000;
  e=1*10^-6;
  k=(ea/dx^2)*[2,-1,0,0,0;-1,2,-1,0,0;0,-1,2,-1,0;0,0,-1,2,-1;0,0,dx/2,-2*dx,3*dx/2];
  A=inv(M)*k;

  [p,v,error,it]= potencia (A,e);
  wn=p^0.5;
  disp('El menor autovalor es'); disp(p);
  disp('El vector asociado es'); disp(v);
  disp('El error es');disp(error);
  disp('El número de iteraciones es'); disp(it);
  disp('w numérica es'); disp(wn);
  
  [r,error,it]=raiz (e);
  disp('El valor de lambda es'); disp(r);
  disp('El error es');disp(error);
  disp('El número de iteraciones es'); disp(it);
  wa=(r/l)*(ea/ra)^0.5;
  disp('w analítica es'); disp(wa);
endfunction

function [p,v,error,it]=potencia (A,ep)
 
y0=[1;1;1;1;1]; 
alfa0=[2;3;1;1;1];
fin=0;
it=0;
Ainv=inv(A);

 while fin==0
  it=it+1;
  x=y0/max(abs(y0));
  y=Ainv*x;
  alfa=y./x;
  error=max(abs(alfa-alfa0)./abs(alfa));
  if error<ep
    fin=1;
  endif
  if it==1000
 fin=2;
endif

y0=y;
alfa0=alfa;
endwhile
  v=y/max(abs(y));
  p=1/max(alfa);
endfunction


function [r,error,it]=raiz (e)

a=0;
b=1;
it=0;
numsol=0;
while numsol==0
  it=it+1;
  r=(a+b)/2;
  error=(abs(fc(r)));
  if error < e
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

function fx= fc (lam)
  fx= lam*tan(lam)-((8*10^-9)*(5*10^4))/1.02;
 endfunction