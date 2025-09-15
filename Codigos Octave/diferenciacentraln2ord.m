function diferenciacentraln2ord
  yA=pi/4;
  tA=0;
  dt=0.01;
  tf=1;
  dyA=0;
  L=10;
  g=9.81;
  yV=yA-dt*dyA+0.5*(dt^2)*((-g/L)*yA);#Calculo de y viejo inicial (serie de Taylor)
  it=0;
  
  for i=tA:dt:tf
    it=it+1;
    yN=(-g/L)*yA*(dt^2)-yV+2*yA;#Ecuacion de Recurrencia (se calcula)
    vx(it)=i;
    vy(it)=yA;
    yV=yA;
    yA=yN;
  endfor
  plot(vx,vy);
endfunction