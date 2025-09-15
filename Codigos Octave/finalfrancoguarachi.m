function repaso
  ri=0; dr=100; rf=600;
  r=(ri+dr):dr:(rf);
  T=25000;
  p=80;
  L=600;
  Aa=100;

  K1=[-2/3, 2/3, 0, 0, 0;
      1, -2, 1, 0, 0;
      0, 1, -2, 1, 0;
      0, 0, 1, -2, 1;
      0, 0, 0, 1, -2];
  K1=(-T/(dr^2))*K1;
  K2=[-4/(3*r(1)), 4/(3*r(1)), 0, 0, 0;
      -1/r(2), 0, 1/r(2), 0, 0;
      0, -1/r(3), 0, 1/r(3), 0;
      0, 0, -1/r(4), 0, 1/r(4);
      0, 0, 0, -1/r(5), 0];
  K2=(-T/(2*dr))*K2;
  A=zeros(5,1);
  A(5)=(-T/(dr^2))*Aa+(-T/(2*dr))*Aa;

  it=0;
  for i=r(2):dr:r(6)
    it=it+1;

    y(it)=g(i);
  endfor

  K=K1+K2;
  f=y'-A;

  #Expresar K y f
  disp(K);
  disp("");
  disp(f);
  inv(K)*f
  #Encontrar la solucion
  D=diag(diag(K));
  B=K-D;
  T=-inv(D)*B;
  c=inv(D)*f;

  u0=ones(5,1);
  ea=0.0000001;
  it2=0;
  fin=0;

  while fin==0
    it2=it2+1;

    uN=T*u0+c;
    err=max(abs((uN-u0)./uN));
    if err<=ea
      fin=1;
    endif
    if it2==1000
      fin=1;
    endif
    u0=uN;
  endwhile
  disp("");
  u=zeros(7,1);
  u(1)=(4/3)*uN(1)-(1/3)*uN(2);
  for i=1:5
    u(i+1)=uN(i);
  endfor
  u(7)=0;
  disp(u);

  #Derivada de u
  du=derivadaprimera(u,dr);
  disp("");
  disp(du);

  #Integral
  r2=ri:dr:rf;
  vec=25000*2*pi*r2*du*du;
  I=TC(vec,dr);
  disp("");
  disp(I);

endfunction
#Funcion
function y=g(x)
  y=80*cos((pi*x)/(2*600));
endfunction
#Derivada Primera
function VD=derivadaprimera(vec,dx)
  n=length(vec);
  VD=zeros(n,1);
  VD(1)=((-3/(2*dx))*vec(1))+(2/dx)*vec(2)-(1/(2*dx))*vec(3);
  for i=2:(n-1)
    VD(i)=(vec(i+1)-vec(i-1))/(2*dx);
  endfor
  VD(n)=(3/(2*dx))*vec(n)-(2/dx)*vec(n-1)+(1/(2*dx))*vec(n-2);
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
