function clase7
A=1;
k=2;
w=5;
xi=0; dx=(10/50); xf=10;
x=xi:dx:xf;
ti=0; dt=(2/100); tf=2;
t=ti:dt:tf;
Y=zeros(51,101);

for i=1:length(x)
  for j=1:length(t)
    Y(i,j)=A*sin(k*x(i)-w*t(j));
  endfor
endfor

#Punto 3
I3=TC(Y(6,51:101),dt);
disp("");
disp("Punto3");
disp(I3);

#Punto 4
I1=TC(Y(11,:),dt);
I2=TC(Y(1:6,16),dx);
I4=I1+I2;
disp("");
disp("Punto4");
disp(I4);

#Punto 5
I5=zeros(10,1);
for i=1:length(x)
  I5(i)=TC(Y(i,:),dt);
endfor
disp("");
disp("Punto5");
disp(sum(I5));

#Punto 6
dy6=derivadaprimera(Y(31,:),dt);
disp("");
disp("Punto6");
disp(sum(dy6));

#Punto 7
dy7=derivadaprimera(Y(:,31),dx);
disp("");
disp("Punto7");
disp(sum(dy7));

#Punto 8
dy1=derivadaprimera(Y(11,:),dx);
I1=TC(dy1,dt);

dy2=derivadaprimera(Y(:,16),dx);
I2=TC(dy2,dx);

dy3=derivadaprimera(Y(:,16),dt);
I3=TC(dy3,dx);

I=I1+I2+I3;
disp("");
disp("Punth8");
disp(I);

#PLOT
Y1=Y(1,:);
Y2=Y(26,:);
Y3=Y(51,:);
hold on
plot(t,Y1,"r");
plot(t,Y2,"g");
plot(t,Y3,"b");
hold off

endfunction
#Integral
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
#Derivada Primera
function VD=derivadaprimera(vy,paso)
  m=length(vy);
  VD(1)=(-3/(2*paso))*vy(1)+(2/(paso))*vy(2)+(-1/(2*paso))*vy(3);
  for i=2:m-1
    VD(i)=(vy(i+1)-vy(i-1))/(2*paso);
  endfor
  VD(m)=(3/(2*paso))*vy(m)+(-2/(paso))*vy(m-1)+(1/(2*paso))*vy(m-2);
endfunction
