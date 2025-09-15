function diferenciacentraln
  yA=2;
  xA=0;
  dx=0.01;
  xf=1;
  yV=yA-dx*(2*yA-2*xA-1);#Calculo de y viejo inicial (serie de Taylor)
  it=0;
  
  for i=xA:dx:xf
    it=it+1;
    yN=(2*dx)*(2*yA-2*i-1)+yV;#Ecuacion de Recurrencia (se calcula)
    vx(it)=i;
    vy(it)=yA;
    yV=yA;
    yA=yN;
  endfor
  plot(vx,vy);
endfunction