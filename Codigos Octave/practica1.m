function repaso
  M=zeros(5,5);
  M(1,1)=1;M(2,2)=4;M(3,3)=8;M(4,4)=12;M(5,5)=16;
  M=5*M;
  K=zeros(5,5);
  for i=1:5
    K(i,i)=2;
  endfor
  K(2,2)=4;
  for i=1:4
    K(i+1,i)=-1;
    K(i,i+1)=-1;
  endfor
  K=5*K;
  
  A=inv(M)*K;
  y0=5*[1;2;3;4;5];
  a0=[1;1;1;1;1];
  ea=10*(10^(-9));
  fin=0;
  it=0;
  
  while fin==0
    it=it+1;
    xn=y0/max(abs(y0));
    yn=A*xn;
    a=yn./xn;
    error=max(abs(a-a0)./abs(a));
    
    if error<=ea
      fin=1;
    endif
    av(:,it)=a0;
    y0=yn;
    a0=a;
  endwhile
  avec=y0;
  aval=1/max(a0);
  #disp(aval);
  paso=pi/100;
  x=1:paso:1+pi;
  vy=avec(1)*(aval+sin(x));
  I=SC(vy,paso);
  disp(I);
endfunction
function I=SC(vy,paso)
  m=length(vy);
  w=zeros(m,1);
  w(1)=1/3;
  w(m)=1/3;
  for i=2:2:m-1
    w(i)=4/3;#impar
  endfor
  for i=3:2:m-2
    w(i)=2/3;#par
  endfor
  w=paso*w;
  I=0;
  for i=1:m
    I=I+w(i)*vy(i);
  endfor
endfunction