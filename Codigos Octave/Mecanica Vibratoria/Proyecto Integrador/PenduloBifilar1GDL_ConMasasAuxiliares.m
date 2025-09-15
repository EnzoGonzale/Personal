L=0.265; D=0.6; mb=0.642; I=0.0343; g=9.81;
m=57*10^(-3);

Y_centimetros=0:5:25;
Y=Y_centimetros*(1/100);

T2=zeros(1,length(Y));
Im=zeros(1,length(Y));

for i=1:length(Y)
T2(i)=(16*pi^2*L*(I+2*m*(Y(i)^2)))/((mb+2*m)*g*D^2);
Im(i)=2*m*(Y(i))^2;
endfor

plot(Im,T2);
title("T^2 vs Im", 'fontsize', 18);
set(gca, 'fontsize', 16);

I_equivalente=zeros(1,length(T2));
for i=1:length(T2)
I_equivalente(i)=(((mb+2*m)*g*D^2)/(16*(pi^2)*L))*T2(i);
endfor
figure();
plot(Y_centimetros,I_equivalente);
title("Inercia equivalente vs Y", 'fontsize', 18);
set(gca, 'fontsize', 16);
disp("Momento de Inercia I_equivalente cuando las masas se encuentran en el centro");
I_equivalente(1)
