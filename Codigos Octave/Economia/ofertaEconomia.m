function ofertaEconomia


  #Por Día
  Q=0:1:250;
  QM=0.5:1:249.5
  CF=600000*(1/24)+13306.01*(1/24)+30010.19*(1/24)+83000/24+42000/24+500000/24+2*1000000/24+700000/24+1000000/24+83333*(1/24);
  CV=31000*15*10^(-3)+7549.08*(50/3600)+300+122.77+300.83*500*10^(-3);
  IM=2400*ones(1,length(Q));

  #Por Mes
  #Q=0:1:4500;
  #QM=0.5:1:4499.5;
  #CF=600000+13306.01+30010.19+83000+42000+500000+2*1000000+700000+1000000+83333;
  #CV=(2108000+226472.4+4512.45+552465);
  #IM=10.8*10^6*ones(1,length(Q));

  CT=CF*ones(1,length(Q))+(CV+0.0007*CV*Q).*Q;
  CTMe=CT./Q;
  CVMe=CV./Q;
  for i=1:length(Q)-1
    CM(i)=CT(i+1)-CT(i);
  endfor

  hold on
  plot(QM,CM,'r');
  #plot(Q,CVMe,'b');
  #plot(Q,CTMe,'g');
  #plot(Q,IM,'y');
  hold off



xlabel("Cantidad ofrecida de cafes por dia");
ylabel("Precio de la unidad del Café");
% Obtener el handle del objeto de texto del xlabel
xlabel_handle = get(gca, 'xlabel');

% Establecer el nuevo tamaño de fuente para xlabel
set(xlabel_handle, 'FontSize', 18); % Cambia al tamaño que desees

% Obtener el handle del objeto de texto del ylabel
ylabel_handle = get(gca, 'ylabel');

% Establecer el nuevo tamaño de fuente para ylabel
set(ylabel_handle, 'FontSize', 18); % Cambia al tamaño que desees

% Obtener el handle de los ejes actuales
axes_handle = gca;

% Establecer el nuevo tamaño de fuente para las etiquetas de los ejes (los números)
set(axes_handle, 'FontSize', 16); % Cambia al tamaño que desees
endfunction
