function newton_raphson
  it=0;
  x0=0;
  fin=0;
  itmax=4;
  while fin==0
    it=it+1;
    r=x0-(fc(x0)/df(x0));
    Dr=r-x0;
    Cr=abs(Dr)/r;
    if it==itmax
      fin=1;
    endif
    disp(it);
    disp(r);
    disp(Dr);
    disp(Cr);
    x0=r;
  endwhile
endfunction
function y=fc(x)
  y=144-360*exp(-2*x)+72*exp(-3*x);
endfunction
function dy=df(x)
  dy=720*exp(-2*x)-216*exp(-3*x);
endfunction