function t2
xi=0;
xf=2*pi;
dx=pi/10;
it=0;
tol=0.2;
for i=xi:dx:xf
  it=it+1;
  x(it)=i;
  y(it)=aproxcos(x(it),tol);
  y2(it)=cos(x(it));
  y3(it)=(y2(it)-y(it))*10;
endfor
r=raices1(xi,xf,tol);
hold on
    plot(x,y,"g");
    plot(x,y2,"*b");
    plot(x,y3,"r");
    plot(r,0,"*y");
hold off
    title("Enzo Gonzalez 13369");
    xlabel("x");
    ylabel("y");
endfunction
#
function F=factor1(n)
  F=1;
  for i=1:n
    if i<=n
    F=F*i;
    endif
  endfor
endfunction
#
function out=aproxcos(x,tol)
  e=999;
  out=0;
  outp=0;
  i=0;
  while e>tol
    out=out+(-1)^i*x^(2*i)/factor1(2*i);
    if i>0
      e=abs(out-outp);
    endif
    outp=out;
    i=i+1;
  endwhile
endfunction
#
function r=raices1(a,b,error1)
  it=0;
  numsol=0;
  tol=0.2;
  while numsol==0
  it=it+1;
  r=0.5*(a+b);
  errorr=abs(aproxcos(r,tol));
  if (errorr<=error1)
    numsol=1;
  endif
  if it==1000
     numsol=2;
  endif
  if (aproxcos(a,tol)*aproxcos(r,tol)>0)
     a=r;
  else
     b=r;
   endif
 endwhile
 disp('La raiz es: ');disp(r);
 disp("numero de it: ");disp(it);
endfunction