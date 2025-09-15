function repaso
  ri=0; rf=3.4; dr=(rf/6);
  r=ri:dr:rf;
  ti=0; tf=90000; dt=100;
  t=ti:dt:tf;

  k=60.5; c=343; p=7850; Te=400; Te0=0; td=2000;

  M=zeros(5);
  for i=1:5
    M(i,i)=c*p;
  endfor

  K1=zeros(5);
  for i=1:5
    K1(i,i)=-2;
  endfor
  for i=1:4
    K1(i,i+1)=1;
    K1(i+1,i)=1;
  endfor
  K1(1,1)=-2/3; K1(1,2)=2/3;
  K1=(k/(dr^2))*K1;

  K2=zeros(5);
  for i=1:4
    K2(i+1,i)=-1/r(i+1);
    K2(i,i+1)=1/r(i);
  endfor
  K2(1,1)=-4/(3*r(2));
  K2(1,2)=4/(3*r(2));
  K2=(k/(2*dr))*K2;

  b1=zeros(5,1);
  b1(5)=(k*Te)/(dr^2)-Te/(2*dr*r(6));

  K=inv(M)*(K1+K2);
  b=inv(M)*b1;

  #Funcion g
  it=0;
  for i=ti:dt:tf
    it=it+1;
    y(it)=g(i);
  endfor

  #Runge Kutta
  it2=0;
  w=0.5;
  u0=zeros(5,1);
  for i=ti:dt:tf
    it2=it2+1;

    du=K*u0+b*g(i);
    k1=dt*du;

    tg=i+(dt/(2*w));
    ug=u0+(k1/(2*w));
    dug=K*ug+b*g(tg);
    k2=dt*dug;

    uN=u0+(1-w)*k1+w*k2;
    m(:,it2)=uN;

    u0=uN;
  endfor
  #Fin de RK
  nr=length(r);
  nt=length(t);
  it3=0;
  u1=m(1,:);
  u2=m(4,:);
  for i=ti:dt:tf
    it3=it3+1;
    u3(it3)=Te*g(i);
  endfor

  hold on
  plot(t,u1,"r");
  plot(t,u2,"b");
  plot(t,u3,"g");
  hold off

  u=zeros(nr,nt);
  u(1,:)=(4/3)*m(1,:)-(1/3)*m(2,:);
  for i=1:5
    u(i+1,:)=m(i,:);
  endfor
  u(nr,:)=u3;

  #Derivada Primera
  du=derivadaprimera(u(:,900),dr);

  #Integral
  vec=2*pi*r.*(du.^2);
  I=TC(vec,dr);
  I
endfunction
#Funcion g
function y=g(x)
  if x<=2000
    y=(1/2000)*x;
  else
    y=1;
  endif
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
#Derivada Primera
function VD=derivadaprimera(vy,paso)
  m=length(vy);
  VD(1)=(-3/(2*paso))*vy(1)+(2/(paso))*vy(2)+(-1/(2*paso))*vy(3);
  for i=2:m-1
    VD(i)=(vy(i+1)-vy(i-1))/(2*paso);
  endfor
  VD(m)=(3/(2*paso))*vy(m)+(-2/(paso))*vy(m-1)+(1/(2*paso))*vy(m-2);
endfunction
