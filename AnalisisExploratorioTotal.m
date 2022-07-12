clear
clc

% Lectura de ficheros
delimeters = {'tab','\t'};
T = readtable('DatosCompletos.txt','Delimiter',delimeters);
T08 = readtable('Datos08.txt','Delimiter',delimeters);
T12 = readtable('Datos12.txt','Delimiter',delimeters);
T44 = readtable('DatosCompletos_BDS44.txt','Delimiter',delimeters);


%% Datos por Avión
icaos = unique(T.ICAO);
icaos_44 = unique(T44.ICAO);
datoAvion = table;
j = 1;


Icao = 322;

for i=1:(length(T.ICAO))
    if ismember(icaos(Icao),T.ICAO(i))
        datoAvion(j,:) = T(i,:);
        j = j + 1;
    end
    if mod(i,50000) == 0
        i
    end
end

%Pasar de TimeStamp a HH:MM:SS en CET
Tiempo = [];
for i=1:length(datoAvion.TimeStamp)
    UTC_epoch_seconds=datoAvion.TimeStamp(i)+3600;
    UTC_offset=UTC_epoch_seconds/(24*60*60);
    atomTime=UTC_offset+datenum(2020,2,23);
    d = datetime(atomTime,'ConvertFrom','datenum');
    Tiempo = [Tiempo;d];
end
fprintf('Fin \n');
datosAvion = [datoAvion,array2table(Tiempo)];

%% Representacion de rutas en el mapa de España
figure
ax = worldmap('spain');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
geoshow(T08.Lat ,T08.Lon,'DisplayType','point','Marker','.');

%% clasificación de datos por altura_08
sop=1:1:15;
Id1=zeros(size(T08,1),length(sop));

for i=1:(length(T08.Alt))
    a = floor(T08.Alt(i)/1000) +1;
    if a <= 15
        Id1(i,a) = 1;
    else 
        Id1(i,:) = 0;
    end
end

%% Boxplot Viento, Temperatura y Presion por alturas_08

Id = fliplr(Id1); %usar solo una vez

figure
subplot(2,2,1)
boxplot(T08.WindSpeed,Id,'orientation','horizontal','Labels',1:15)
title('Velocidad del viento por altura')
xlabel('Velocidad del viento [m/s]')
ylabel('Altura [km]')
xlim([0 60])
grid on

subplot(2,2,2)
boxplot(T08.WindDirec_,Id,'orientation','horizontal','Labels',1:15)
title('Dirección del viento por altura')
xlabel('Dirección del viento [°]')
ylabel('Altura [km]')
grid on

subplot(2,2,3)
boxplot(T08.Temperature,Id,'orientation','horizontal','Labels',1:15)
title('Temperatura por altura')
xlabel('Temperatura [K]')
ylabel('Altura [km]')
grid on

subplot(2,2,4)
boxplot(T08.Pressure,Id,'orientation','horizontal','Labels',1:15)
title('Presión por altura')
xlabel('Presión [Pa]')
ylabel('Altura [km]')
ax = gca;
ax.XAxis.Exponent = 0;
grid on

%% Numero de datos por altura_08
datosCol = zeros(1,15);
for i = 1:15
    Ceros = find(~Id(:,i));
    datosCol(1,i) = length(T08.Alt)-length(Ceros);
end
% stairs(datosCol)

histogram(T.Alt/1000,0:1:15)
xlim([0 16])
xticks([0:16])
title('Número de datos por altura')
ylabel('Número de datos')
xlabel('Altura [km]')
grid on

%% Clasificar datos por horas

sot=1:1:9;
h=zeros(size(T,1),length(sot));

for i=1:(length(T.TimeStamp))
    a = floor(T.TimeStamp(i)/3600) - 7;
    h(i,a) = 1;
end

datosHora = zeros(1,9);
for i = 1:9
    Ceros = find(~h(:,i));
    datosHora(1,i) = length(T.TimeStamp)-length(Ceros);
end

t1 = datetime(2020,2,23,9,30,0);
t2 = datetime(2020,2,23,17,30,0);
t = t1:hours(1):t2;

bar(t,datosHora,1)
title('Número de datos por hora')
ylabel('Número de datos')
xlabel('Horas')
grid on



%% Altura y Temperatura de un Avión con el tiempo

%Altura
figure

subplot(1,2,1)
plot(datosAvion.Tiempo,datosAvion.Alt)
title('Altura del avión a lo largo del tiempo:',icaos(Icao))
ylabel('Altura [km]')
xlabel('Hora')
grid on

%Temperatura
subplot(1,2,2)
plot(datosAvion.Tiempo,datosAvion.Temperature)
title('Temperatura que mide el avión a lo largo del tiempo:',icaos(Icao))
ylabel('Temperatura [K]')
xlabel('Hora')
grid on

%% Altura, Temperatura y Roll de un Avión con el tiempo

%Altura
figure
subplot(2,2,[1,3])
plot(datosAvion.Tiempo,datosAvion.Alt)
title('Altura del avión a lo largo del tiempo:',icaos(Icao))
ylabel('Altura [km]')
xlabel('Hora')
grid on

%Temperatura
subplot(2,2,2)
plot(datosAvion.Tiempo,datosAvion.Temperature)
title('Temperatura que mide el avión a lo largo del tiempo:',icaos(Icao))
ylabel('Temperatura [K]')
xlabel('Hora')
grid on

%Roll
subplot(2,2,4)
plot(datosAvion.Tiempo,datosAvion.Roll)
title('Roll del avión a lo largo del tiempo:',icaos(Icao))
ylabel('Roll [º]')
xlabel('Hora')
grid on

%% Trayectoria de un avión en 3D sobre mapa
uif = uifigure;
g = geoglobe(uif);
geoplot3(g,datosAvion.Lat,datosAvion.Lon,datosAvion.Alt,'y','LineWidth',2) 

%%
geoshow(datosAvion.Lat,datosAvion.Lon)
axis equal
%% Posición de un avión con velocidad del viendo y heading:

u = -0.02 .* datosAvion.WindSpeed .* sind(datosAvion.WindDirec_);
v = -0.02 .* datosAvion.WindSpeed .* cosd(datosAvion.WindDirec_);
u_ = 1 .* sind(datosAvion.TrackAngle);
v_ = 1 .* cosd(datosAvion.TrackAngle);

figure 
hold on
for i=1:length(u)
    if mod(i,9) == 0
        quiver(datosAvion.Lon(i),datosAvion.Lat(i), u(i),v(i),'r')
        hold on
        quiver(datosAvion.Lon(i),datosAvion.Lat(i), u_(i),v_(i),'b')
    end
end
hold off
% 
% quiver(datosAvion.Lon,datosAvion.Lat, u,v)
% hold on
% quiver(datosAvion.Lon,datosAvion.Lat, u_,v_)
% grid on
legend ('wind', 'direction')
title('Trayectoria 2D del avión:',icaos(Icao))
ylabel('Latitud [º]')
xlabel('Longitud [º]')
axis equal
hold off


%% Rosa de los vientos por altura

polarbubblechart(datosAvion.WindDirec_,datosAvion.Alt,datosAvion.WindSpeed)

%%
% scatter3(datosAvion.WindSpeed,datosAvion.Alt,datosAvion.Mach)
% 
% xlabel('Velocidad Viento [m/s]')
% ylabel('Altitud [m]')
% zlabel('Mach')


%% quitar outliers, no funciona con tablas pero si con vectores
Temp = rmoutliers(Temp,'quartiles');
WS = rmoutliers(WS,'quartiles');
Press = rmoutliers(Press,'quartiles');

%% Aviones con datos y clasificación de datos por avión

datoAvion44 = table;
j = 1;


Icao_44 = 4; % Del 1 al 7

for i=1:(length(T44.ICAO))
    if ismember(icaos_44(Icao_44),T44.ICAO(i))
        datoAvion44(j,:) = T44(i,:);
        j = j + 1;
    end
end

%Pasar de TimeStamp a HH:MM:SS en CET
Tiempo44 = [];
for i=1:length(datoAvion44.TimeStamp)
    UTC_epoch_seconds=datoAvion44.TimeStamp(i)+3600;
    UTC_offset=UTC_epoch_seconds/(24*60*60);
    atomTime=UTC_offset+datenum(2020,2,23);
    d = datetime(atomTime,'ConvertFrom','datenum');
    Tiempo44 = [Tiempo44;d];
end
fprintf('Fin \n');
datosAvion44 = [datoAvion44,array2table(Tiempo44)];

%%
figure
plot(datosAvion44.Tiempo44,datosAvion44.Temperature)
title('Temperatura que mide el avión a lo largo del tiempo:',icaos_44(Icao_44))
ylabel('Temperatura [K]')
xlabel('Hora')
grid on

figure
plot(datosAvion.Tiempo,datosAvion.Temperature)
title('Temperatura derivada del avión a lo largo del tiempo:',icaos(Icao))
ylabel('Temperatura [K]')
xlabel('Hora')
grid on

%%
