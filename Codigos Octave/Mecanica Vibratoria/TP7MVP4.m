#Punto 4 TP7 M. Vibratoria
k=1000; m=100; z=0.1; Tf=2*pi;
w=(2*pi)/Tf;
wn=sqrt(k/m);
wd=wn*sqrt(1-z^2);
B=w/wn;
dt=(2/190)*pi;
t=0:dt:3*pi;

a1=-4/(pi^2); a2=0; a3=-4/(9*pi^2); a4=0; a5=-4/(25*pi^2);

b1=0; b2=0; b3=0; b4=0; b5=0;

x=zeros(1,length(t));
P=zeros(1,length(t));
PF=zeros(1,length(t));

for i=1:5
  H(i)=(1/k)*(1/sqrt((1-(i*B)^2)^2+(2*z*i*B)^2));
  if (i*B)>1
  phi(i)=pi-atan((2*z*i*B)/(1-(i*B)^2));
  else
  phi(i)=atan((2*z*i*B)/(1-(i*B)^2));
  endif

endfor
a5=2*a5; a4=2*a4;
a1=2*a1; a2=2*a2; a3=2*a3;

for i=1:length(t)

  t_mod = mod(t(i),Tf);# Calcula la posición dentro de un solo periodo (entre 0 y 2*pi)
  if (t_mod>=0 && t_mod<=pi)
    P(i)=(4/(2*pi))*(t_mod)-1;
  endif
  if (t_mod>pi && t_mod<=2*pi)
    P(i)=1-(4/(2*pi))*(t_mod-((2*pi)/2));
  endif


  x(i)=a1*H(1)*cos(1*w*t(i)-phi(1))+a2*H(2)*cos(2*w*t(i)-phi(2))+a3*H(3)*cos(3*w*t(i)-phi(3))+a4*H(4)*cos(4*w*t(i)-phi(4))+a5*H(5)*cos(5*w*t(i)-phi(5))+b1*H(1)*sin(1*w*t(i)-phi(1))+b2*H(2)*sin(2*w*t(i)-phi(2))+b3*H(3)*sin(3*w*t(i)-phi(3))+b4*H(4)*sin(4*w*t(i)-phi(4))+b5*H(5)*sin(5*w*t(i)-phi(5));
  PF(i)=a1*cos(w*t(i))+a2*cos(2*w*t(i))+a3*cos(3*w*t(i))+a4*cos(4*w*t(i))+a5*cos(5*w*t(i))+b1*sin(1*w*t(i))+b2*sin(2*w*t(i))+b3*sin(3*w*t(i))+b4*sin(4*w*t(i))+b5*sin(5*w*t(i));

endfor

plot(t,x);
figure();
hold on
plot(t,PF);
plot(t,P);
hold off
