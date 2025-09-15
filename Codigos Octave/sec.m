function sec
  it=0;
  itmax=1000;
  x1=0; x2=1;
  ea=0.00001;
  numsol=0;
  while numsol==0
  it=it+1;
  m=(fcc(x2)-fcc(x1))/(x2-x1);
  xr=x2-(fcc(x2)/m);
  errorr=abs((x2-x1)/x2);
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
 disp('La raiz es: ');disp(xr);
 disp("numero de it: ");disp(it);
 disp("El error es: ");disp(errorr);
 
 hold on
 plot(e,"r");
 plot(xrr,"g");
 hold off
endfunction
#
function fun=fcc(x)
  fun=exp(x)+2^(-x)-6+2*cos(x);
endfunction