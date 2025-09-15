function repaso
  ri=0; dr=(500/7); rf=500;
  r=ri:dr:rf;
  T=25000; m=0.8;

  M=zeros(6);
  for i=1:6
    M(i,i)=m;
  endfor

  K1=zeros(6);
  for i=1:6
    K1(i,i)=-2;
  endfor
  for i=1:5
    K1(i+1,i)=1;
    K1(i,i+1)=1;
  endfor
  K1(1,1)=-2/3; K1(1,2)=2/3;
  K1=(-T/(dr^2))*K1;

  K2=zeros(6);
  for i=1:5
    K2(i+1,i)=-1/r(i+2);
    K2(i,i+1)=1/r(i+1);
  endfor
  K2(1,1)=-4/(3*r(2)); K2(1,2)=4/(3*r(2));
  K2=(-T/(2*dr))*K2;

  K=K1+K2;

  A=inv(M)*K;

  #Potencia
  a0=ones(6,1);
  y0=ones(6,1);
  it=0;
  fin=0;
  ea=0.0001;

  while fin==0
    it=it+1;

    xn=y0/max(abs(y0));
    yn=inv(A)*xn;
    an=yn./xn;

    err=max(abs((an-a0)./an));
    if err<=ea
      fin=1;
    endif
    if it==1000
      fin=1;
    endif
    a0=an;
    y0=yn;
  endwhile
  aval=1/max(an)
  avec=yn/max(yn)
  u=zeros(8,1);
  u(1)=(4/3)*yn(1)-(1/3)*yn(2);
  for i=1:6
    u(i+1)=yn(i);
  endfor

  #Aproximacion
  Fi=zeros(6,2);
  for i=1:6
    Fi(i,1)=f1(r(i+1));
    Fi(i,2)=f2(r(i+1));
  endfor
  a=(Fi'*Fi)*(Fi'*yn);

  it2=0;
  for i=r(2):dr:r(7)
    it2=it2+1;

    g(it2)=a(1)*f1(i)+a(2)*f2(i);
  endfor
  r1=r(2):dr:r(7);

  g=g/max(g);
  figure(1);
  plot(r,u);
  figure(2);
  plot(r1,g);

  #Integral
  du=derivadaprimera(u,dr)
  vec=du'*T*du*2*pi*r;
  I=TC(vec,dr)
endfunction
function y1=f1(x)
  y1=cos((pi*x)/(2*500));
endfunction
function y2=f2(x)
  y2=cos((pi*x)/(4*500));
endfunction
function I=TC(vec,dx)
  n=length(vec);
  m=ones(n,1);
  m(1)=0.5;
  m(n)=0.5;
  m=dx*m;
  I=0;
  for i=1:n
    I=I+m(i)*vec(i);
  endfor
endfunction
function VD=derivadaprimera(vec,dx)
  n=length(vec);
  VD=zeros(n,1);
  VD(1)=(vec(2)-vec(1))/dx;
  for i=2:n-1
    VD(i)=(vec(i+1)-vec(i-1))/(2*dx);
  endfor
  VD(n)=(vec(n)-vec(n-1))/dx;
endfunction
