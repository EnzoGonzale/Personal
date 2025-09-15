function graficas_economia
  Q=0:1:250;
  GramosDeCafe=0:32:8000;
  CF=(50000+6267.745)*ones(1,length(Q));
  CostoDeCadaCafe = (5700/30)+24.5*(24033/500)+3050;
  CV = ones(1,length(Q))*1.5*CostoDeCadaCafe;

  CT = CF + CV.*Q;
  #CM=CT./Q;
  CM=zeros(1,length(Q));
  CM(1)=CT(1);
  for i=1:length(Q)-1
    CM(i+1) = CT(i+1)-CT(i);
  endfor
  plot(Q,CM);
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
