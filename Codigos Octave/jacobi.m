function jacobi
  A=[2.5,-1,0;-1,3,-2;0,-1,3.5];
  b=[1;2;3];
  y0=[1;1;1];
  ep=0.000001;
  it=0;
  itmax=1000;
  fin=0;
  
  D=diag(diag(A));
  B=A-D;
  T=-inv(D)*B;
  c=inv(D)*b;
  
  while fin==0
    it=it+1;
    yn=T*y0+c;
    error=max(abs(yn-y0)./abs(yn));
    if(error<ep)
    fin=1;
  endif
  if(it==itmax)
  fin=2;
endif
my(:,it)=y0;
y0=yn;
endwhile
hold on
plot(my(1,:));
plot(my(2,:));
plot(my(3,:));
hold off
disp("error "); disp(error);
disp("iteraciones "); disp(it);
disp("vector "); disp(yn);
endfunction
