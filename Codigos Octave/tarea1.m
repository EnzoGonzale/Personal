function tarea1
  disp("Legajo:"); n=input("");
  if 9000<n && n<14000
    it=0;
    for i=0:9
      it=it+1;
      x(it)=i;
      y(it)=fc(x(it),n);
    endfor
    hold on
    plot(x,y,"*r");
    hold off
    title("FUNCIONES");
    xlabel("x");
    ylabel("y");
  else
    disp("El numero de Legajo no corresponde a un alumno de CNYC 2022.");
  endif
endfunction
function y=fc(x,N)
  y=x*x+2*x+(1/N);
endfunction