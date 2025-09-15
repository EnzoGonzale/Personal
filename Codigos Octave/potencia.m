function potencia
  A=[2.5,-1,0;-1,3,-2;0,-1,3.5];
  y0=[1;1;1];
  a0=[2;3;1];#alfa inicial
  ep=0.000001;
  it=0;
  itmax=1000;
  fin=0;
  
  Ai=inv(A);
  
  while fin==0
    it=it+1;
    xn=y0/max(abs(y0));
    yn=Ai*xn;
    a=yn./xn;#alfa nuevo
    error=max(abs(a-a0)./abs(a));
    if(error<ep)
    fin=1;
  endif
  if(it==itmax)
  fin=2;
endif
av2(:,it)=a0;
y0=yn;
a0=a;
v=yn;
nv=max(abs(v));
av=1/max(a);
endwhile
hold on
plot(av2');
hold off
disp("norma del auto vector ");disp(nv);
disp("menor autovalor "); disp(av);
endfunction