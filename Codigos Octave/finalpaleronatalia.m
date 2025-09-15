function repaso
  ri=0; dr=(90/6); rf=90;
  r=ri:dr:rf;

  T=490; p=100;

  K1=zeros(5);
  K2=zeros(5);
  for i=1:5
    K1(i,i)=-2;
  endfor
  for i=1:4
    K1(i+1,i)=1;
    K1(i,i+1)=1;

    K2(i+1,i)=-1/r(i+2);
    K2(i,i+1)=1/r(i+1);
  endfor
  K1(1,1)=-2/3; K1(1,2)=2/3;
  K1=(T/(dr^2))*K1;
  K2(1,1)=-4/(3*r(2)); K2(1,2)=4/(3*r(2));
  K2=(T/(2*dr))*K2;
  K=K1+K2

  it=0;
  for i=r(2):dr:r(6)
    it=it+1;
    q(it)=f(i);
  endfor
  q

  #Jacobi
  D=diag(diag(K));
  B=K-D;
  TT=-inv(D)*B;
  C=inv(D)*q';
  fin=0;
  u0=zeros(5,1);
  ea=0.0001;

  while fin==0
    un=TT*u0+C;
    err=max(abs((un-u0)./un));
    if err<=ea
      fin=1;
    endif
    u0=un;
  endwhile

  u=zeros(7,1);
  u(1)=(4/3)*un(1)-(1/3)*un(2);
  for i=1:5
    u(i+1)=un(i);
  endfor
  u

  #Derivada Primera
  du=derivadaprimera(u,dr);
  figure(1);
  plot(r,du);

  #Inregral
  vI=du'*T*du*2*pi*r;
  I=TC(vI,dr);
  I=0.5*I

  #Runge Kutta
  y0=5;
  B=3;
  L=10; g=1000;
  w=0.5;
  it2=0;

  for i=0:0.5:1
    it2=it2+1;

    dy=B;
    dB=-(g/L)*sin(y0);
    k1=0.5*dy;
    k2=0.5*dB;

    tg=i+(0.5/(2*w));
    yg=y0+(k1/(2*w));
    Bg=B+(k2/(2*w));
    dyg=Bg;
    dBg=-(g/L)*sin(yg);
    k1g=0.5*dyg;
    k2g=0.5*dBg;

    yN=y0+(1-w)*k1+w*k1g;
    BN=B+(1-w)*k2+w*k2g;
    m(it2,1)=yN;

    y0=yN;
    B=BN;
  endfor
  dm=derivadaprimera(m,0.5);
  figure(2);
  plot(m,dm);
  #Fin de RK
endfunction
function y=f(x)
  y=100*cos((pi*x)/(2*90));
endfunction
#Derivada Primera
function VD=derivadaprimera(vec,dx)
  n=length(vec);
  VD=zeros(n,1);
  VD(1)=((-3/2)*vec(1)+2*vec(2)-(1/2)*vec(3))/dx;
  for i=2:n-1
    VD(i)=(vec(i+1)-vec(i-1))/(2*dx);
  endfor
  VD(n)=((3/2)*vec(n)-2*vec(n-1)+(1/2)*vec(n-2))/dx;
endfunction
#Trapecio Compuesto
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
