#Punto 1 Enzo Gonzalez 13369 Industrial
function parcial21
  A=zeros(4,4);
  A(1,1)=52;
  A(2,2)=104;
  A(3,3)=104;
  A(4,4)=52;
  for i=1:3
    A(i+1,i)=-26;
    A(i,i+1)=-26;
  endfor
  b=[-52;104;104;-52];
  y0=zeros(4,1);
  D=diag(diag(A));
  B=A-D;
  T=-inv(D)*B;
  c=inv(D)*b;
  
  it=0;
  itmax=4;
  fin=0;
  
  while fin==0
    it=it+1;
    yn=T*y0+c;
    Dx=yn-y0;
    Cr=max(abs(Dx))/max(abs(yn));
    V=Dx.*yn;
    if it==itmax
      fin=1;
    endif
    disp("Numero de iteracion ");disp(it);
    disp("vector ");disp(yn);
    disp("Dx ");disp(Dx);
    disp("Cr ");disp(Cr);
    disp("vector x solucion ");disp(V);
    disp("-----------------------");
    y0=yn;
  endwhile
endfunction

#Punto2
function parcial22
  it=0;
  itmax=3;
  fin=0;
  x0=0;
  
  while fin==0
    it=it+1;
    r=x0-(fc(x0)/df(x0));
    Draiz=r-x0;
    Cr=abs(Draiz)/abs(r);
    if it==itmax
      fin=1;
    endif
    disp("iteracion");disp(it);
    disp("Raiz ");disp(r);
    disp("Draiz ");disp(Draiz);
    disp("Cr ");disp(Cr);
    disp("------------------");
    x0=r;
  endwhile
endfunction
function y=fc(x)
  y=150-375*exp(-2*x)+75*exp(-3*x);
endfunction
function dy=df(x)
  dy=750*exp(-2*x)-225*exp(-3*x);
endfunction