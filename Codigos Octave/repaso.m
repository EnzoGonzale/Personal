function repaso
  xi=0; dx=(50000/5); xf=50000;
  x=xi:dx:xf;
  EA=24*10^4; pA=8*10^(-9); Ml=1.02*10^(-2); L=50000;

  #Metodo Directo
  Fi=zeros(3);
  x1=0:(L/2):L;
  y1=[0;1;1];
  for i=1:3
    Fi(i,1)=1;
    Fi(i,2)=x1(i);
    Fi(i,3)=(x1(i))^2;
  endfor
  a1=inv(Fi)*y1;
  #Fin MD

  M=zeros(5);
  K1=zeros(5);
  K2=zeros(5);
  for i=1:4
    M(i,i)=pA*f(x(i+1),a1);
    K1(i+1,i)=1;
    K1(i,i+1)=1;
    K1(i,i)=-2;
  endfor
  M(5,5)=-Ml;
  K1(5,4)=0;
  K1=(-EA/(dx^2))*K1;
  K2(5,3)=1/2; K2(5,4)=-2; K2(5,5)=3/2;
  K2=(-EA/dx)*K2;

  K=K1+K2
  M

  #Potencia Inversa
  A=inv(M)*K;
  y0=ones(5,1);
  a0=ones(5,1);
  ea=0.001;
  it=0;
  fin=0;

  while fin==0
    it=it+1;

    xn=y0/max(abs(y0));
    yn=(inv(A))*xn;
    an=yn./xn;
    err=max(abs(((an-a0)./(an))));

    if err<=ea
      fin=1;
    endif
    if it=1000
      fin=1;
    endif
    m(:,it)=yn;
    y0=yn;
    a0=an;
  endwhile
  v=yn/yn(1)
  landa=1/max(an/an(1));
  w=landa^0.5

  u=zeros(6,1);
  u(2:6)=v

  #Derivada
  du=deriprimera(u,dx)

  #Integral
  it2=0;
  for i=xi:dx:xf
    it2=it2+1;

    fy(it2)=f(i,a1);
  endfor
  vec=du*fy*du;
  I=TC(vec,dx);
  En=0.5*I

endfunction
function y=f(x,a)
  y=a(1)*1+a(2)*x+a(3)*(x^2);
endfunction
#Derivada
function VD=deriprimera(v,dx)
  n=length(v);
  VD=zeros(n,1);
  VD(1)=((-3/2)*v(1)+2*v(2)-(1/2)*v(3))*(1/dx);
  for i=2:n-1
    VD(i)=((v(i+1)-v(i-1))/(2*dx));
  endfor
  VD(n)=((3/2)*v(n)-2*v(n-1)+(1/2)*v(n-2))*(1/dx);
endfunction
#Integral
function I=TC(v,dx)
  n=length(v);
  m=ones(n,1);
  m(1)=0.5;
  m(n)=0.5;
  m=dx*m;
  I=0;
  for i=1:n
    I=I+m(i)*v(i);
  endfor
endfunction
