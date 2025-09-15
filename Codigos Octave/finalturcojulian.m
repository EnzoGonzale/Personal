function repaso
  ri=0; dr=(100/2); rf=100;
  r=ri:dr:rf;
  ti=0; dt=0.002; tf=6;
  t=ti:dt:tf;
  T=15000/13; m=0.1; p=10; ro=5;

  w2=((T*2)/((dr^2)*3*m))+((T*4)/(r(2)*dr*2*3*m));
  f=p/m;

  #Runge Kutta
  it=0;
  u0=0;
  B=0;
  w=0.5;

  for i=ti:dt:tf
    it=it+1;

    du=B;
    dB=-w2*u0+f*sin(ro*i);
    k1=dt*du;
    k2=dt*dB;

    tg=i+(dt/(2*w));
    ug=u0+(k1/(2*w));
    Bg=B+(k2/(2*w));
    dug=Bg;
    dBg=-w2*ug+f*sin(ro*tg);
    k1g=dt*dug;
    k2g=dt*dBg;

    uN=u0+(1-w)*k1+w*k1g;
    BN=B+(1-w)*k2+w*k2g;

    u(1,it)=uN;

    u0=uN;
    B=BN;
  endfor
  #Fin de RK

  #Graficar u(L/2,t) y cuantas veces pasa por cero
  figure(1);
  plot(t,u);
  it2=0;
  cc=0;
  for i=ti:dt:tf-1
    it2=it2+1;
    if u(it2)*u(it2+1)<=0;;
      cc=cc+1;
    endif
  endfor
  disp("Veces que se anula la funcion");disp(cc);

  #Derivada Primera
  dut=derivadaprimera(u,dt);
  figure(2);
  plot(u,dut);

  #Integral
  it3=0;
  for i=ti:dt:tf
    it3=it3+1;
    s(it3)=f*sin(ro*i);
  endfor
  vec=s-w2*u;
  I=TC(vec,dt);
  disp("El valor de la integral"); disp(I);

  #
  n=length(t);
  j1=u(n)
  j2=(4/3)*u(n)
  j3=0;
  it4=0;

  for i=ri:dr:rf
    it4=it4+1;
    uu(it4)=g(i);
  endfor
  figure(3);
  plot(r,uu);

  uuu=[j2;j1;0];
  figure(4);
  plot(r,uuu);
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
#Interpolacion
function y=g(x)
  y=17.393*(((x-50)*(x-100))/((17.393-50)*(17.393-100)))+13.045*(((x-0)*(x-100))/((13.045-0)*(13.045-100)));
endfunction
