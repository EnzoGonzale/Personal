function pol1
  vx=[-3;0;5;10];
  vy=[-82.5;48;-94.5;288];
  x=-3:0.1:10;
  for i=1:length(x)
    y(i)=poli(x(i),vx,vy);
    y2(i)=metododirecto(x(i),vx,vy);
  endfor
  hold on
  plot(x,y);
  plot(x,y2);
  hold off
endfunction
#aproximacion
function p=poli(x,vx,vy)
    mf(:,1)=vx.^(3);
    mf(:,2)=vx;
    mf(:,3)=1;
  coef=inv(mf'*mf)*(mf'*vy);
  p=0;
  for i=1:3
    p=p+coef(i)*x^(i-1);
  endfor
endfunction
#interpolacion
function p2=metododirecto(x,vx,vy)
  n=length(vx);
  for i=1:n
    mf(:,i)=vx.^(i-1);
  endfor
  coef=inv(mf)*vy;
  p2=0;
  for i=1:n
    p2=p2+coef(i)*x^(i-1);
  endfor
endfunction