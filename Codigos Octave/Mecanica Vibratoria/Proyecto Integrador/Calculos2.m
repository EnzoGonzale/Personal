#Variables
m1=0.642; m2=0.353; mt=m1+m2;
L1=0.265; L2=0.15; L_barra_grande=0.795; L_barra_chica=0.6; #L1=0.265; L2=0.15;
D1=0.6; D2=0.21; #D1=0.6; D2=0.21;
I1=0.0343; I2=7.125*10^(-3);
g=9.81;
dt=0.01; t=0:dt:40;
F=zeros(1,length(t));

Meq=zeros(2,2); Meq(1,1)=I1; Meq(2,2)=I2;
Keq=zeros(2,2); Keq(1,2)=Keq(2,1)=-(m2/L2)*D1*D2; Keq(1,1)=((mt/L1)+(m2/L2))*D1^2; Keq(2,2)=(m2/L2)*D2^2;
Keq=(g/4)*Keq;

u=I2/I1;

[Xn, lamdas]=eig(Keq,Meq);
w=lamdas.^0.5;

D2=sqrt((4*I2*L2)/(m2*g))*(w(2,2)/(1+u));


u
a
D2
m2
