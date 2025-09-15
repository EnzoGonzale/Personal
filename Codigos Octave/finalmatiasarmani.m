function repaso
  ti=0; dt=1/10; tf=150; td=50;
  t=ti:dt:tf;
  it=0;
  #Funcion g
  for i=ti:dt:tf
    it=it+1;

    y(it)=g(i);
  endfor

  K=zeros(5);
  for i=1:5
    K(i,i)=2;
  endfor
  for i=1:4
    K(i,i+1)=-1;
    K(i+1,i)=-1;
  endfor

  b=[5;8;9;8;5];
  b=(50/18)*b;

  M=zeros(5);
  for i=1:5
    M(i,i)=1;
  endfor

  #Diferencia Central
  zA=zeros(5,1);
  dz=zeros(5,1);
  zV=zeros(5,1);
  it2=0;

  for i=ti:dt:tf
    it2=it2+1;

    zN=2*zA-zV+(dt^2)*(b*g(i)-K*zA);
    m(:,it2)=zN;
    zV=zA;
    zA=zN;
  endfor

  #Graficar z3
  z3=m(3,:);
  figure(1);
  plot(t,z3);

  #Calcular y Graficar v3
  v3=derivadaprimera(z3,dt);
  figure(2);
  plot(t,v3);

  #Aproximacion
  t2=td:dt:tf;
  z32=z3(501:1501); %Usando m=(t/dt)+1
  n2=length(t2);
  Fi=ones(n2,3);
  it3=0;
  it4=0;
  for i=td:dt:tf
    it3=it3+1;
    Fi(it3,2)=f1(i);
    Fi(it3,3)=f2(i);
  endfor
  a=inv(Fi'*Fi)*Fi'*z32';
  for i=td:dt:tf
    it4=it4+1;

    Fa(it4)=a(1)*1+a(2)*f1(i)+a(3)*f2(i);
  endfor
  disp("vector Fi'*Fi"); disp(Fi'*Fi);
  disp("vector de coeficientes"); disp(a);
  figure(3);
  plot(t2,Fa);

  dif=abs(z32-Fa);
  figure(4);
  plot(t2,dif);
endfunction
#Funcion g
function y=g(x)
  if x<=50
    y=(1/50)*x;
  else
    y=1;
  endif
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
#Funciones Aproximacion
function y1=f1(x)
  y1=sin((2*pi*x)/27.142);
endfunction
function y2=f2(x)
  y2=cos((2*pi*x)/27.142);
endfunction
