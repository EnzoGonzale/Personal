function rungekutta
  x0=0;
  y0=2;
  dx=0.01;
  xf=1;
  w=0.5;
  it=0;
  x=x0:dx:xf;
  for i=x0:dx:xf
    it=it+1;
    
    dy=2*y0-2*i-1;
    k1=dx*dy;
    
    xg=i+(dx/(2*w));
    yg=y0+(k1/(2*w));
    dyg=2*yg-2*xg-1;
    k2=dx*dyg;
    
    vx(it)=i;
    vy(it)=y0;
    
    y0=y0+(1-w)*k1+w*k2;
  endfor
  plot(vx,vy);
  endfunction