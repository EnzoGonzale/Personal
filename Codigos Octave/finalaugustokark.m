function repaso
  xi=0; dx=(0.55/5); xf=0.55;
  x=xi:dx:xf;
  ti=0; dt=10; tf=5000;
  t=ti:dt:tf;

  T=22.798*10^(-5); p=300; td=1000;

  K=zeros(4);
  for i=1:4
    K(i,i)=-2;
  endfor
  for i=1:3
    K(i+1,i)=1;
    K(i,i+1)=1;
  endfor
  K=(T/(dx^2))*K

  b=zeros(4,1);
  b(1)=(T/(dx^2))*p

  #Runge Kutta
  it=0;
  w=0.5;
  u0=ones(4,1);

  for i=ti:dt:tf
    it=it+1;

    du=K*u0+b*g(i);
    k1=dt*du;

    tg=i+(dt/(2*w));
    ug=u0+(k1/(2*w));
    dug=K*ug+b*g(tg);
    k2=dt*dug;

    uN=u0+(1-w)*k1+w*k2;
    m(:,it)=uN;

    u0=uN;
  endfor
  #Fin de RK

  #Graficos de ui y uL
  nt=length(t);
  nx=length(x);
  it2=0;
  u=zeros(nx,nt);
  for i=ti:dt:tf
    it2=it2+1;
    u(1,it2)=p*g(i);
  endfor
  for i=1:4
    u(i+1,:)=m(i,:);
  endfor
  u(6,:)=(4/3)*m(4,:)-(1/3)*m(3,:);

  ui=u(1,:);
  uL=u(nx,:);
  hold on
  plot(t,ui,"r");
  plot(t,uL,"b");
  hold off

  #Expresar u en t=3000
  ut=u(:,301)

  #Calcular la derivada en t=3000
  du=derivadaprimera(ut,dx)
  figure(2);
  plot(x,du);

  #Integral
  vec=du.^2;
  I=TC(vec,dx)
endfunction
function y=g(x)
  if x<=1000
    y=(1/1000)*x;
  else
    y=1;
  endif
endfunction
#Derivada Primera
function VD=derivadaprimera(vec,dx)
  n=length(vec);
  VD(1)=((-3/2)*vec(1)+2*vec(2)+(1/2)*vec(3))*(1/dx);
  for i=2:n-1
    VD(i)=(vec(i+1)-vec(i-1))/(2*dx);
  endfor
  VD(n)=((3/2)*vec(n)-2*vec(n-1)+(1/2)*vec(n-2))*(1/dx);
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
