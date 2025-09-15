function repaso
  T=10000;
  dr=50/3;
  r0=0; L=100;
  r=dr:dr:L-dr;

  K1=(-T/(dr^2))*[-2/3, 2/3, 0, 0, 0;
                  1, -2, 1, 0, 0;
                  0, 1, -2, 1, 0;
                  0, 0, 1, -2, 1;
                  0, 0, 0, 1, -2];
  K2=(-T/(2*dr))*[-4/(3*r(1)), 4/(3*r(1)), 0,0,0;
                  -1/(r(2)), 0, 1/(r(2)), 0, 0;
                  0, -1/r(3), 0, 1/r(3), 0;
                  0, 0, -1/r(4), 0, 1/r(4);
                  0, 0, 0, -1/r(5), 0];
  K=K1+K2;
  M=zeros(5,5);
  for i=1:5
    M(i,i)=0.8;
  endfor
  A=inv(M)*K;

  fin=0;
  it=0;
  ea=0.00001;
  y0=ones(5,1);
  a0=ones(5,1);

  while fin==0
    it=it+1;

    xn=y0/max(abs(y0));
    yn=inv(A)*xn;
    an=yn./xn;
    err=max(abs((an-a0)./an));

    if it==1000
      fin=1;
    endif
    if err<=ea
      fin=1;
    endif

    avec=yn/max(abs(yn));
    aval=1/max(an);
    aval=(aval^(1/2));

    y0=yn;
    a0=an;
  endwhile

  disp(avec);
  disp(aval);
  disp(it);

  u=zeros(7,1);
  u(1)=(4/3)*avec(1)-(1/3)*avec(2);
  for i=1:5
    u(i+1)=avec(i);
  endfor
  disp(u);

  #Derivada
  VD=derivadaprimera(u,dr);
  VD=VD';
  disp(VD);

  #Integral
  rr=r0:dr:L;
  vec=T*2*pi*rr*VD*VD;
  I=TC(vec,dr);
  disp(I);

  #Aprox
  it2=0;
  it3=0;
  for i=r0:dr:L
    it2=it2+1;
    y1(it2,1)=f1(i);
  endfor
  for i=r0:dr:L
    it3=it3+1;
    y2(it3,1)=f2(i);
  endfor
  Fi=zeros(7,2);
  for i=1:7
    Fi(i,1)=y1(i);
    Fi(i,2)=y2(i);
  endfor
  veca=inv(Fi'*Fi)*Fi'*u;
  polaprox=veca(1)*f1(rr)+veca(2)*f2(rr);
  hold on
  plot(rr,polaprox, "r");
  plot(rr,u, "g*");
  hold off
endfunction
#Derivada
function VD=derivadaprimera(vec,dx)
  m=length(vec);
  VD(1)=((-3/2)*vec(1)+(2*vec(2))+(-1/2)*vec(3))/dx;
  for i=2:m-1
    VD(i)=(vec(i+1)-vec(i-1))/(2*dx);
  endfor
  VD(m)=((3/2)*vec(m)+(-2*vec(m-1))+((1/2)*vec(m-2)))/dx;
endfunction
#Integral
function I=TC(vec,dx)
  m=length(vec);
  w=ones(m,1);
  w(1)=0.5;
  w(m)=0.5;
  for i=1:m
    I=I+dx*vec(i)*w(i);
  endfor
endfunction
function y1=f1(x)
  y1=cos((pi*x)/200);
endfunction
function y2=f2(x)
  y2=cos((pi*3*x)/200);
endfunction
