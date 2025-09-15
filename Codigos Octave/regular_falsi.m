function regular_falsi
  it=0;
  itmax=8;
  xi=0;
  xf=5;
  fin=true;
  while fin==true
    it=it+1;
    m=(fc(xf)-fc(xi))/(xf-xi);
    r=xf-(fc(xf)/m);
    if it==itmax
      fin=false;
    endif
    disp(it);
    disp(r);
    if fc(r)*fc(xf)>0
      xf=r;
  else
    xi=r;
    endif
  endwhile
endfunction
function y=fc(x)
  y=3*(x^2)-12;
endfunction