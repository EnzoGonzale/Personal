mb=0.642; m=0.050; mt=mb+m;
L1=0.265; L2=0.15; L_barra=0.795; #L1=0.265; L2=0.15;
D=0.6;#D1=0.6;
I=0.0343;
g=9.81;
dt=0.01; t=0:dt:40;
F=zeros(1,length(t));
k1=(mt*g)/(2*L1);
k2=(m*g)/L2;

M=zeros(2); M(1,1)=I+(m*D^2)/2; M(2,2)=m*(L2^2); #M(2,1)=M(1,2)=-m*D*L2;
C=zeros(length(M(1,:))); #C(2,1)=C(1,2)=-m*D*L2;
K=zeros(length(M(1,:))); K(2,1)=K(1,2)=-k2*D*L2; K(1,1)=((D^2)/2)*(k1+k2); K(2,2)=2*k2*(L2^2);

w_barra=sqrt(K(1,1)/M(1,1));

u=M(2,2)/M(1,1);
u
u_optimo=0.1;
u_optimo
m_optima=(u*M(1,1))/2;
m_optima

a=1/(1+u);
wTMD=w_barra*a;

L2_optimo=(2*g)/(wTMD^2);
a
L2_optimo


u2=0.01;
m2=u2*mb*0.5;
a2=1/(1+u2);
wTMD2=w_barra*a2;
L22=(2*g)/(wTMD2^2);

u2
m2
a2
wTMD2
L22


