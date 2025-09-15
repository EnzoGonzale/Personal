function gaussseidel
   A=zeros(4,4);
  for i=1:3
    A(i+1,i)=-2.5;
    A(i,i+1)=-2.5;
  endfor
  for i=1:4
    A(i,i)=10;
  endfor
  x0=[0;0;0;0];
  b=[25;25;25;25];
  D=diag(diag(A));
  B=A-D;
  T=-inv(D)*B;
  TS=triu(T);
  TI=tril(T);
  c=inv(D)*b;
  it=0;
  itmax=10;
  fin=true;
  while fin==true;
    it=it+1;
    yn=TI*x0+TS*x0+c;
    if it==itmax
      fin=false;
    endif
    avec=max(abs(yn));
    disp(yn);
    disp("\n");
    x0=yn;
  endwhile
endfunction