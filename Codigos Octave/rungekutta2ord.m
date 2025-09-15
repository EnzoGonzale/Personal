function rungekutta2ord
  t0=0;
  dt=0.01;
  tf=1;
  y1=pi/4;
  y2=0;
  L=10;
  g=9.81;
  w=0.5;
  it=0;
  
  for i=t0:dt:tf
    it=it+1;
    dy1=y2;
    dy2=-(g/L)*y1;
    k1=dt*dy1;
    k2=dt*dy2;
    
    tg=i+(dt/(2*w));
    y1g=y1+(k1/(2*w));
    y2g=y2+(k2/(2*w));
    dy1g=y2g;
    dy2g=-(g/L)*dy1g;
    k1g=dt*dy1g;
    k2g=dt*dy2g;
    
    vt(it)=i;
    vy1(it)=y1;
    
    y1=y1+(1-w)*k1+w*k1g;
    y2=y2+(1-w)*k2+w*k2g;
  endfor
  plot(vt,vy1);
endfunction