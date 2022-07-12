clear
clc

delimeters = {'tab','\t'};
datos = readtable('DatosCompletos.txt','Delimiter',delimeters); %%Datos completos derivados
era5 = 'ERA5_data.nc';

datosTotal = readtable('TotalDatos_Avion_Reanalisis.txt','Delimiter',delimeters); %%Datos derivados + reanálisis


time = ncread(era5,'time');
time = transpose((time -1053168)*3600); %tiempo en formato TimeStamp

tempK = ncread(era5,'t');
tempC = tempK - 273.15;

longitud = ncread(era5,'longitude');
latitud = ncread(era5,'latitude');

press = ncread(era5,'level'); %27 niveles de presión en milibares
press = press * 100; %27 niveles de presión en Pa

u = ncread(era5,'u');
v = ncread(era5,'v');

%% Temperatura por altura para primera hora
t0_1= transpose(tempC(:,:,27,1));
t1_1= transpose(tempC(:,:,24,1));
t2_1= transpose(tempC(:,:,21,1));
t3_1= transpose(tempC(:,:,17,1));
t4_1= transpose(tempC(:,:,15,1));
t5_1= transpose(tempC(:,:,14,1));
t6_1= transpose(tempC(:,:,12,1));
t7_1= transpose(tempC(:,:,11,1));
t8_1= transpose(tempC(:,:,10,1));
t9_1= transpose(tempC(:,:,9,1));
t10_1= transpose(tempC(:,:,8,1));
t11_1= transpose(tempC(:,:,8,1));
t12_1= transpose(tempC(:,:,7,1));
t13_1= transpose(tempC(:,:,6,1));
t14_1= transpose(tempC(:,:,5,1));
t15_1= transpose(tempC(:,:,4,1));


%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VISUALIZACION ERA5 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gráficas de Temperatura por altura para primera hora
load coastlines.mat

figure
colormap turbo(70)


subplot(4,4,1)
mymap = pcolor(longitud,latitud,t0_1)
mymap.EdgeAlpha = 0.05;
mymap.CDataMapping = 'scaled';
title('0 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar('eastoutside')
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,2)
mymap = pcolor(longitud,latitud,t1_1)
mymap.EdgeAlpha = 0.05;
title('1 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,3)
mymap = pcolor(longitud,latitud,t2_1)
mymap.EdgeAlpha = 0.05;
title('2 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,4)
mymap = pcolor(longitud,latitud,t3_1)
mymap.EdgeAlpha = 0.05;
title('3 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,5)
mymap = pcolor(longitud,latitud,t4_1)
mymap.EdgeAlpha = 0.05;
title('4 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,6)
mymap = pcolor(longitud,latitud,t5_1)
mymap.EdgeAlpha = 0.05;
title('5 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,7)
mymap = pcolor(longitud,latitud,t6_1)
mymap.EdgeAlpha = 0.05;
title('6 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,8)
mymap = pcolor(longitud,latitud,t7_1)
mymap.EdgeAlpha = 0.05;
title('7 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,9)
mymap = pcolor(longitud,latitud,t8_1)
mymap.EdgeAlpha = 0.05;
title('8 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,10)
mymap = pcolor(longitud,latitud,t9_1)
mymap.EdgeAlpha = 0.05;
title('9 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,11)
mymap = pcolor(longitud,latitud,t10_1)
mymap.EdgeAlpha = 0.05;
title('10 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,12)
mymap = pcolor(longitud,latitud,t11_1)
mymap.EdgeAlpha = 0.05;
title('11 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,13)
mymap = pcolor(longitud,latitud,t12_1)
mymap.EdgeAlpha = 0.05;
title('12 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,14)
mymap = pcolor(longitud,latitud,t13_1)
mymap.EdgeAlpha = 0.05;
title('13 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,15)
mymap = pcolor(longitud,latitud,t14_1)
mymap.EdgeAlpha = 0.05;
title('14 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

subplot(4,4,16)
mymap = pcolor(longitud,latitud,t15_1)
mymap.EdgeAlpha = 0.05;
title('15 km')
xlabel('Longitud [º]')
ylabel('Latitud [º]')
hold on
plot(coastlon,coastlat,'k')
colorbar
hold off
ax = gca;
ax.CLim = [-75 35];

%% Temperatura por altura y por hora
load coastlines.mat
for i=1:10
    t0_1= transpose(tempC(:,:,27,i));
    t1_1= transpose(tempC(:,:,24,i));
    t2_1= transpose(tempC(:,:,21,i));
    t3_1= transpose(tempC(:,:,17,i));
    t4_1= transpose(tempC(:,:,15,i));
    t5_1= transpose(tempC(:,:,14,i));
    t6_1= transpose(tempC(:,:,12,i));
    t7_1= transpose(tempC(:,:,11,i));
    t8_1= transpose(tempC(:,:,10,i));
    t9_1= transpose(tempC(:,:,9,i));
    t10_1= transpose(tempC(:,:,8,i));
    t11_1= transpose(tempC(:,:,8,i));
    t12_1= transpose(tempC(:,:,7,i));
    t13_1= transpose(tempC(:,:,6,i));
    t14_1= transpose(tempC(:,:,5,i));
    t15_1= transpose(tempC(:,:,4,i));
    
    
    
    figure
    colormap turbo(90)
    
    subplot(4,4,1)
    mymap = pcolor(longitud,latitud,t0_1)
    mymap.EdgeAlpha = 0.05;
    mymap.CDataMapping = 'scaled';
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,2)
    mymap = pcolor(longitud,latitud,t1_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,3)
    mymap = pcolor(longitud,latitud,t2_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,4)
    mymap = pcolor(longitud,latitud,t3_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,5)
    mymap = pcolor(longitud,latitud,t4_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,6)
    mymap = pcolor(longitud,latitud,t5_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,7)
    mymap = pcolor(longitud,latitud,t6_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,8)
    mymap = pcolor(longitud,latitud,t7_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,9)
    mymap = pcolor(longitud,latitud,t8_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,10)
    mymap = pcolor(longitud,latitud,t9_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,11)
    mymap = pcolor(longitud,latitud,t10_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,12)
    mymap = pcolor(longitud,latitud,t11_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,13)
    mymap = pcolor(longitud,latitud,t12_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,14)
    mymap = pcolor(longitud,latitud,t13_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,15)
    mymap = pcolor(longitud,latitud,t14_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(4,4,16)
    mymap = pcolor(longitud,latitud,t15_1)
    mymap.EdgeAlpha = 0.05;
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];

end


    %% Temperatura en una única altura por horas  
load coastlines.mat

    t0_1= transpose(tempC(:,:,24,2));
    t1_1= transpose(tempC(:,:,24,3));
    t2_1= transpose(tempC(:,:,24,4));
    t3_1= transpose(tempC(:,:,24,5));
    t4_1= transpose(tempC(:,:,24,6));
    t5_1= transpose(tempC(:,:,24,7));
    t6_1= transpose(tempC(:,:,24,8));
    t7_1= transpose(tempC(:,:,24,9));
    t8_1= transpose(tempC(:,:,24,10));
    figure
    colormap turbo(90)
    
    subplot(3,3,1)
    mymap = pcolor(longitud,latitud,t0_1)
    mymap.EdgeAlpha = 0.05;
    title('09:00')
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(3,3,2)
    mymap = pcolor(longitud,latitud,t1_1)
    mymap.EdgeAlpha = 0.05;
    title('10:00')
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(3,3,3)
    mymap = pcolor(longitud,latitud,t2_1)
    mymap.EdgeAlpha = 0.05;
    title('11:00')
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(3,3,4)
    mymap = pcolor(longitud,latitud,t3_1)
    mymap.EdgeAlpha = 0.05;
    title('12:00')
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(3,3,5)
    mymap = pcolor(longitud,latitud,t4_1)
    mymap.EdgeAlpha = 0.05;
    title('13:00')
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(3,3,6)
    mymap = pcolor(longitud,latitud,t5_1)
    mymap.EdgeAlpha = 0.05;
    title('14:00')
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(3,3,7)
    mymap = pcolor(longitud,latitud,t6_1)
    mymap.EdgeAlpha = 0.05;
    title('15:00')
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(3,3,8)
    mymap = pcolor(longitud,latitud,t7_1)
    mymap.EdgeAlpha = 0.05;
    title('16:00')
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];
    
    subplot(3,3,9)
    mymap = pcolor(longitud,latitud,t8_1)
    mymap.EdgeAlpha = 0.05;
    title('17:00')
    hold on
    plot(coastlon,coastlat,'k')
    colorbar
    hold off
    ax = gca;
    ax.CLim = [-75 35];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% COMPARACION DAT vs ERA5 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
T_era5 = zeros(1,length(datos.ICAO));
u_era5 = zeros(1,length(datos.ICAO));
v_era5 = zeros(1,length(datos.ICAO));
for i=1:(length(datos.ICAO))
    %Comparamos el instante temporal del mensaje:
    VecDifTiempo = zeros(1,length(time));
    for j=1:length(time)
        VecDifTiempo(j) = abs(datos.TimeStamp(i)-time(j));
    end
    [difTiempo,PosTime] = min(VecDifTiempo);
    
    %Comparamos longitud:
    VecDifLon = zeros(1,length(longitud));
    for j=1:length(longitud)
        if datos.Lon(i)>180
            datos.Lon(i) = datos.Lon(i) - 360;
            if longitud(j)>0
                VecDifLon(j) = abs(datos.Lon(i)-longitud(j));
            else
                VecDifLon(j) =  abs(datos.Lon(i)+longitud(j));
            end
        else
            if longitud(j)>0
                VecDifLon(j) = abs(datos.Lon(i)+longitud(j));
            else
                VecDifLon(j) = abs(datos.Lon(i)-longitud(j));
            end
        end
    end
    [difLon,PosLon] = min(VecDifLon);
    LonAlrededor = zeros(1,3);
    if PosLon ~= 1 || PosLon ~= 51
        r=1;
        for j=PosLon-1:PosLon+1
            LonAlrededor(r) = abs(datos.Lon(i)-longitud(j));
            r=r+1;
        end
        [DifLonLevel, LonCercano] = min([abs(LonAlrededor(1)-LonAlrededor(2)) abs(LonAlrededor(2)-LonAlrededor(3))]);
        if LonCercano == 1
            LonCercano = -1;
        else
            LonCercano = 1;
        end
    elseif PosLon == 1
        LonCercano = 1;
    elseif PosLon == 51
        LonCercano = -1;
    end

    %Comparamos latitud:
    VecDifLat = zeros(1,length(latitud));
    for j=1:length(latitud) %no hace falta diferencia signo por que españa está siempre con latitud +
        VecDifLat(j) = abs(datos.Lat(i)-latitud(j));
    end
    [difLat,PosLat] = min(VecDifLat);
    LatAlrededor = zeros(1,3);
    if PosLat ~= 1 || PosLat ~= 41
        r=1;
        for j=PosLat-1:PosLat+1
            LatAlrededor(r) = abs(datos.Lat(i)-latitud(j));
            r=r+1;
        end
        [DifLatLevel, LatCercano] = min([abs(LatAlrededor(1)-LatAlrededor(2)) abs(LatAlrededor(2)-LatAlrededor(3))]);
        if LatCercano == 1
            LatCercano = -1;
        else
            LatCercano = 1;
        end
    elseif PosLon == 1
        LonCercano = 1;
    elseif PosLon == 41
        LonCercano = -1;
    end


    %Comparamos el nivel del mensaje (presión):
    VecDifPress = zeros(1,length(press));
    for j=1:length(press)
        VecDifPress(j) = abs(datos.Pressure(i)-press(j));
    end
    [difPress,PosPress] = min(VecDifPress);

    PressAlrededor = zeros(1,3);
    if PosPress ~= 1 || PosPress ~= 10
        r=1;
        for j=PosPress-1:PosPress+1
            PressAlrededor(j) = abs(datos.Pressure(i)-press(j));
            r=r+1;
        end
        [DifPressLevel, NivelCercano] = min([abs(PressAlrededor(1)-PressAlrededor(2)) abs(PressAlrededor(2)-PressAlrededor(3))]);
        if NivelCercano == 1
            NivelCercano = -1;
        else
            NivelCercano = 1;
        end
    elseif PosLon == 1
        LonCercano = 1;
    elseif PosLon == 10
        LonCercano = -1;
    end

    %Vamos a hacer las interpolaciones, en total 8 por variable-> a:"nivel actual", c:"nivel cercano"

    P_a = press(PosPress);
    P_c = press(PosPress+NivelCercano);

    T_1a = tempK(PosLon,PosLat,PosPress,PosTime);
    T_2a = tempK(PosLon,PosLat+LatCercano,PosPress,PosTime);
    T_3a = tempK(PosLon+LonCercano,PosLat,PosPress,PosTime);
    T_4a = tempK(PosLon+LonCercano,PosLat+LatCercano,PosPress,PosTime);
    T_1c = tempK(PosLon,PosLat,PosPress+NivelCercano,PosTime);
    T_2c = tempK(PosLon,PosLat+LatCercano,PosPress+NivelCercano,PosTime);
    T_3c = tempK(PosLon+LonCercano,PosLat,PosPress+NivelCercano,PosTime);
    T_4c = tempK(PosLon+LonCercano,PosLat+LatCercano,PosPress+NivelCercano,PosTime);

    u_1a = u(PosLon,PosLat,PosPress,PosTime);
    u_2a = u(PosLon,PosLat+LatCercano,PosPress,PosTime);
    u_3a = u(PosLon+LonCercano,PosLat,PosPress,PosTime);
    u_4a = u(PosLon+LonCercano,PosLat+LatCercano,PosPress,PosTime);
    u_1c = u(PosLon,PosLat,PosPress+NivelCercano,PosTime);
    u_2c = u(PosLon,PosLat+LatCercano,PosPress+NivelCercano,PosTime);
    u_3c = u(PosLon+LonCercano,PosLat,PosPress+NivelCercano,PosTime);
    u_4c = u(PosLon+LonCercano,PosLat+LatCercano,PosPress+NivelCercano,PosTime);

    v_1a = v(PosLon,PosLat,PosPress,PosTime);
    v_2a = v(PosLon,PosLat+LatCercano,PosPress,PosTime);
    v_3a = v(PosLon+LonCercano,PosLat,PosPress,PosTime);
    v_4a = v(PosLon+LonCercano,PosLat+LatCercano,PosPress,PosTime);
    v_1c = v(PosLon,PosLat,PosPress+NivelCercano,PosTime);
    v_2c = v(PosLon,PosLat+LatCercano,PosPress+NivelCercano,PosTime);
    v_3c = v(PosLon+LonCercano,PosLat,PosPress+NivelCercano,PosTime);
    v_4c = v(PosLon+LonCercano,PosLat+LatCercano,PosPress+NivelCercano,PosTime);

    %Hacemos una primera interpolación en latitud:
    latDifTotal = abs(latitud(PosLat+LatCercano)-latitud(PosLat));
    latDifPeque = abs(datos.Lat(i)-latitud(PosLat));
    y = latDifPeque/latDifTotal;

    T_11a = T_1a*(1-y)+T_2a*y;
    T_21a = T_3a*(1-y)+T_4a*y;
    T_11c = T_1c*(1-y)+T_2c*y;
    T_21c = T_3c*(1-y)+T_4c*y;

    u_11a = u_1a*(1-y)+u_2a*y;
    u_21a = u_3a*(1-y)+u_4a*y;
    u_11c = u_1c*(1-y)+u_2c*y;
    u_21c = u_3c*(1-y)+u_4c*y;

    v_11a = v_1a*(1-y)+v_2a*y;
    v_21a = v_3a*(1-y)+v_4a*y;
    v_11c = v_1c*(1-y)+v_2c*y;
    v_21c = v_3c*(1-y)+v_4c*y;

    %Hacemos una segunda interpolación en longitud:
    lonDifTotal = abs(longitud(PosLat+LonCercano)-longitud(PosLat));
    lonDifPeque = abs(datos.Lon(i)-longitud(PosLat));
    x = lonDifPeque/lonDifTotal;
    
    T_12a = T_11a*(1-x)+T_21a*x;
    T_12c = T_11c*(1-x)+T_21c*x;

    u_12a = u_11a*(1-x)+u_21a*x;
    u_12c = u_11c*(1-x)+u_21c*x;

    v_12a = v_11a*(1-x)+v_21a*x;
    v_12c = v_11c*(1-x)+v_21c*x;

    %Hacemos la última iteración entre niveles de presión:
    pressDifTotal = abs(double(P_a)-double(P_c));
    pressDifPeque = abs(datos.Pressure(i)-double(P_a));
    z = pressDifPeque/pressDifTotal;
    
    T_era5(i) = T_12a*(1-z)+T_12c*z;

    u_era5(i) = u_12a*(1-z)+u_12c*z;

    v_era5(i) = v_12a*(1-z)+v_12c*z;
    
    if mod(i,5000) == 0
        i
    end

end

%Guardamos los datos de era5 al final de la tabla:
    datosTotal = datos;
    datosTotal.Temp_era5 = transpose(T_era5);
    datosTotal.u_era5 = transpose(u_era5);
    datosTotal.v_era5 = transpose(v_era5);

%%
datosTotal.u = (datosTotal.WindSpeed .* sind(datosTotal.WindDirec_));
datosTotal.v = (datosTotal.WindSpeed .* cosd(datosTotal.WindDirec_));

datosTotal = movevars(datosTotal,'u','After','WindDirec_');
datosTotal = movevars(datosTotal,'v','After','u');
%%
writetable(datos,'TotalDatos_Avion_Reanalisis.txt','Delimiter','\t','WriteRowNames',true);
type TotalDatos_Avion_Reanalisis.txt;

%% Comprobación de precisión de datos
dif_Temp = zeros(1,length(datosTotal.ICAO));
dif_u = zeros(1,length(datosTotal.ICAO));
dif_v = zeros(1,length(datosTotal.ICAO));
for i=1:length(datosTotal.ICAO)
    dif_Temp(i)= abs(datosTotal.Temperature(i)-datosTotal.Temp_era5(i));
    dif_u(i)= abs((datosTotal.WindSpeed(i) .* sind(datosTotal.WindDirec_(i)))-datosTotal.u_era5(i));
    dif_v(i)= abs((datosTotal.WindSpeed(i) .* cosd(datosTotal.WindDirec_(i)))-datosTotal.v_era5(i));
end

error_medio_temp = mean(dif_Temp)
error_medio_u = mean(dif_u)
error_medio_v = mean(dif_v)

%Juntar datos
datosTotal.error_Temp = datosTotal.Temperature-datosTotal.Temp_era5;
datosTotal.error_u = (datosTotal.WindSpeed .* sind(datosTotal.WindDirec_))-datosTotal.u_era5;
datosTotal.error_v = (datosTotal.WindSpeed .* cosd(datosTotal.WindDirec_))-datosTotal.v_era5;

%%
writetable(datosTotal,'TotalDatos_Avion_Reanalisis2.txt','Delimiter','\t','WriteRowNames',true);
type TotalDatos_Avion_Reanalisis2.txt;

%% Clasificación datos por altura
sop=1:1:15;
Id=zeros(size(datosTotal,1),length(sop));

for i=1:(length(datosTotal.Alt))
    a = floor(datosTotal.Alt(i)/1000) + 1;
    if a <= 15
        Id(i,a) = 1;
    else 
        Id(i,:) = 0;
    end
end

%% Boxplot Viento, Temperatura y Presion por alturas_08

Id2 = fliplr(Id);
 
figure
subplot(1,3,1)
boxplot(datosTotal.error_Temp,Id2,'orientation','horizontal','Labels',1:15)
title('Error medio de la temperatura por alturas')
xlabel('Grados [º]')
ylabel('Altura [km]')
yticklabels(0:1:14)
grid on

subplot(1,3,2)
boxplot(datosTotal.error_u,Id2,'orientation','horizontal','Labels',1:15)
title('Error medio de la componente u por alturas')
xlabel('Velocidad del viento [m/s]')
ylabel('Altura [km]')
yticklabels(0:1:14)
grid on

subplot(1,3,3)
boxplot(datosTotal.error_v,Id2,'orientation','horizontal','Labels',1:15)
title('Error medio de la componente v por alturas')
xlabel('Velocidad del viento [m/s]')
ylabel('Altura [km]')
yticklabels(0:1:14)
grid on

%% Media de valores por niveles de presión
vecMedia_T_era5 = zeros(10,27); %Temperatura media:10 filas de horas y 27 columnas de presiones
for i=1:10
    for j = 1:27
        vecMedia_T_era5(i,j)=mean(mean(tempK(:,:,j,i)));
    end
end

%Eliminamos la primera y última fila de las que no tendremos info en los datos de navegación
vecMedia_T_era5(1,:)=[];
vecMedia_T_era5(9,:)=[];

MED_T_ERA5 = mean(vecMedia_T_era5);
%% Media de valores por niveles de presión
med_T_press1 = [];
med_T_press2 = [];
med_T_press3 = [];
med_T_press4 = [];
med_T_press5 = [];
med_T_press6 = [];
med_T_press7 = [];
med_T_press8 = [];
med_T_press9 = [];
med_T_press10 = [];
med_T_press11 = [];
med_T_press12 = [];
med_T_press13 = [];
med_T_press14 = [];
med_T_press15 = [];
med_T_press16 = [];
med_T_press17 = [];
med_T_press18 = [];
med_T_press19 = [];
med_T_press20 = [];
med_T_press21 = [];
med_T_press22 = [];
med_T_press23 = [];
med_T_press24 = [];
med_T_press25 = [];
med_T_press26 = [];
med_T_press27 = [];

for i=1:(length(datosTotal.ICAO))
    if (datosTotal.Pressure(i)<=11125)
        med_T_press1 = [med_T_press1 datosTotal.Temperature(i)];
    elseif (11125<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=13250)
        med_T_press2 = [med_T_press2 datosTotal.Temperature(i)];
    elseif (13250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=16250)
        med_T_press3 = [med_T_press3 datosTotal.Temperature(i)];
    elseif (16250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=18750)
        med_T_press4 = [med_T_press4 datosTotal.Temperature(i)];
    elseif (18750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=21250)
        med_T_press5 = [med_T_press5 datosTotal.Temperature(i)];
    elseif (21250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=23750)
        med_T_press6 = [med_T_press6 datosTotal.Temperature(i)];
    elseif (23750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=27500)
        med_T_press7 = [med_T_press7 datosTotal.Temperature(i)];
    elseif (27500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=32500)
        med_T_press8 = [med_T_press8 datosTotal.Temperature(i)];
    elseif (32500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=37500)
        med_T_press9 = [med_T_press9 datosTotal.Temperature(i)];
    elseif (37500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=42500)
        med_T_press10 = [med_T_press10 datosTotal.Temperature(i)];
    elseif (42500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=47500)
        med_T_press11 = [med_T_press11 datosTotal.Temperature(i)];
    elseif (47500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=52500)
        med_T_press12 = [med_T_press12 datosTotal.Temperature(i)];
    elseif (52500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=57500)
        med_T_press13 = [med_T_press13 datosTotal.Temperature(i)];
    elseif (57500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=62500)
        med_T_press14 = [med_T_press14 datosTotal.Temperature(i)];
    elseif (62500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=67500)
        med_T_press15 = [med_T_press15 datosTotal.Temperature(i)];
    elseif (67500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=72500)
        med_T_press16 = [med_T_press16 datosTotal.Temperature(i)];
    elseif (72500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=76250)
        med_T_press17 = [med_T_press17 datosTotal.Temperature(i)];
    elseif (76250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=78750)
        med_T_press18 = [med_T_press18 datosTotal.Temperature(i)];
    elseif (78750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=81250)
        med_T_press19 = [med_T_press19 datosTotal.Temperature(i)];
    elseif (81250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=83750)
        med_T_press20 = [med_T_press20 datosTotal.Temperature(i)];
    elseif (83750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=86250)
        med_T_press21 = [med_T_press21 datosTotal.Temperature(i)];
    elseif (86250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=88750)
        med_T_press22 = [med_T_press22 datosTotal.Temperature(i)];
    elseif (88750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=91250)
        med_T_press23 = [med_T_press23 datosTotal.Temperature(i)];
    elseif (91250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=93750)
        med_T_press24 = [med_T_press24 datosTotal.Temperature(i)];
    elseif (93750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=96250)
        med_T_press25 = [med_T_press25 datosTotal.Temperature(i)];
    elseif (96250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=98750)
        med_T_press26 = [med_T_press26 datosTotal.Temperature(i)];
    elseif (98750<=datosTotal.Pressure(i)) 
        med_T_press27 = [med_T_press27 datosTotal.Temperature(i)];
    end

    if mod(i,5000) == 0
        i
    end
end

MED_T_DATA = [mean(med_T_press1) mean(med_T_press2) mean(med_T_press3) mean(med_T_press4) mean(med_T_press5) mean(med_T_press6) mean(med_T_press7) mean(med_T_press8) mean(med_T_press9) mean(med_T_press10) mean(med_T_press11) mean(med_T_press12) mean(med_T_press13) mean(med_T_press14) mean(med_T_press15) mean(med_T_press16) mean(med_T_press17) mean(med_T_press18) mean(med_T_press19) mean(med_T_press20) mean(med_T_press21) mean(med_T_press22) mean(med_T_press23) mean(med_T_press24) mean(med_T_press25) mean(med_T_press26) mean(med_T_press27)]

%%
error_med_T = [MED_T_DATA;MED_T_ERA5];

error_med_T(:,[1,2,3,4,26,27]) = [];
%%
error_temps = zeros(1,length(error_med_T));
for i = 1:length(error_med_T)
    error_temps(i)=(error_med_T(1,i)-error_med_T(2,i));
end

%% Clasificación datos por niveles de presión:
sop=1:1:27;
Id3=zeros(size(datosTotal,1),length(sop));

for i=1:(length(datosTotal.ICAO))
    if (datosTotal.Pressure(i)<=11125)
        Id3(i,1) = 1;
    elseif (11125<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=13250)
        Id3(i,2) = 1;
    elseif (13250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=16250)
        Id3(i,3) = 1;
    elseif (16250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=18750)
        Id3(i,4) = 1;
    elseif (18750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=21250)
        Id3(i,5) = 1;
    elseif (21250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=23750)
        Id3(i,6) = 1;
    elseif (23750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=27500)
        Id3(i,7) = 1;
    elseif (27500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=32500)
        Id3(i,8) = 1;
    elseif (32500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=37500)
        Id3(i,9) = 1;
    elseif (37500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=42500)
        Id3(i,10) = 1;
    elseif (42500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=47500)
        Id3(i,11) = 1;
    elseif (47500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=52500)
        Id3(i,12) = 1;
    elseif (52500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=57500)
        Id3(i,13) = 1;
    elseif (57500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=62500)
        Id3(i,14) = 1;
    elseif (62500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=67500)
        Id3(i,15) = 1;
    elseif (67500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=72500)
        Id3(i,16) = 1;
    elseif (72500<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=76250)
        Id3(i,17) = 1;
    elseif (76250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=78750)
        Id3(i,18) = 1;
    elseif (78750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=81250)
        Id3(i,19) = 1;
    elseif (81250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=83750)
        Id3(i,20) = 1;
    elseif (83750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=86250)
        Id3(i,21) = 1;
    elseif (86250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=88750)
        Id3(i,22) = 1;
    elseif (88750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=91250)
        Id3(i,23) = 1;
    elseif (91250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=93750)
        Id3(i,24) = 1;
    elseif (93750<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=96250)
        Id3(i,25) = 1;
    elseif (96250<=datosTotal.Pressure(i)) &&  (datosTotal.Pressure(i)<=98750)
        Id3(i,26) = 1;
    elseif (98750<=datosTotal.Pressure(i)) 
        Id3(i,27) = 1;
    else 
        Id3(i,:) = 0;
    end

    if mod(i,10000) == 0
        i
    end
end
%% Media de errores
figure
subplot(1,3,1)
boxplot(datosTotal.error_Temp,Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Error medio de la temperatura por Presión')
xlabel('Grados [º]')
ylabel('Presión [Pa]')
grid on

subplot(1,3,2)
boxplot(datosTotal.error_u,Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Error medio de la componente u por Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

subplot(1,3,3)
boxplot(datosTotal.error_v,Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Error medio de la componente v por Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on
%% Visualización de datos ERA5 en eje z
Id2 = fliplr(Id);

figure
subplot(1,3,1)
boxplot(datosTotal.Temp_era5,Id2,'orientation','horizontal','Labels',1:15)
title('Temperatura por Altura')
xlabel('Temperatura [K]')
ylabel('Altura [km]')
grid on

subplot(1,3,2)
boxplot(datosTotal.u_era5,Id2,'orientation','horizontal','Labels',1:15)
title('Velocidad de la componente u del viento por Altura')
xlabel('Velocidad del viento [m/s]')
ylabel('Altura [km]')
grid on

subplot(1,3,3)
boxplot(datosTotal.v_era5,Id2,'orientation','horizontal','Labels',1:15)
title('Dirección de la componente v del viento por Altura')
xlabel('Velocidad del viento [m/s]')
ylabel('Altura [km]')
grid on


%% Media de datos
figure
subplot(1,3,1)
boxplot(datosTotal.Temperature,Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Temperatura media por nivel de Presión')
xlabel('Grados [º]')
ylabel('Presión [Pa]')
grid on

subplot(1,3,2)
boxplot((datosTotal.WindSpeed .* sind(datosTotal.WindDirec_)),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Componente u media por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

subplot(1,3,3)
boxplot((datosTotal.WindSpeed .* cosd(datosTotal.WindDirec_)),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Componente v media por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on



%% Trabajar con velocidad del viento y dirección en vez de u y v:
WS_era5 = zeros(1,length(datosTotal.ICAO));
WD_era5 = zeros(1,length(datosTotal.ICAO));
for i=1:length(datosTotal.ICAO)
    WS_era5(i) = sqrt((datosTotal.u_era5(i))^2+(datosTotal.v_era5(i))^2);
    if datosTotal.u_era5(i)>0 && datosTotal.v_era5(i)>0
        WD_era5(i)= atand(datosTotal.u_era5(i)/datosTotal.v_era5(i));
    elseif datosTotal.u_era5(i)<0 && datosTotal.v_era5(i)>0
        WD_era5(i)= 360+atand(datosTotal.u_era5(i)/datosTotal.v_era5(i));
    elseif datosTotal.u_era5(i)>0 && datosTotal.v_era5(i)<0
        WD_era5(i)= 180+atand(datosTotal.u_era5(i)/datosTotal.v_era5(i));
    elseif datosTotal.u_era5(i)<0 && datosTotal.v_era5(i)<0
        WD_era5(i)= 180+atand(datosTotal.u_era5(i)/datosTotal.v_era5(i));
    end
end

datosTotal.WS_era5 = transpose(WS_era5);
datosTotal.WD_era5 = transpose(WD_era5);

%% WS y WD ERA5 por presión
figure
subplot(1,2,1)
boxplot((datosTotal.WS_era5),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Velocidad del viento era5 por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

subplot(1,2,2)
boxplot((datosTotal.WD_era5),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Dirección del viento era5 por nivel Presión')
xlabel('Dirección del viento [º]')
ylabel('Presión [Pa]')
grid on

%% WS y WD ERA5 por altura
figure
subplot(1,2,1)
boxplot((datosTotal.WS_era5),Id2,'orientation','horizontal','Labels',0:14)
title('Velocidad del viento ERA5 por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

subplot(1,2,2)
boxplot((datosTotal.WD_era5),Id2,'orientation','horizontal','Labels',0:14)
title('Dirección del viento ERA5 por nivel Presión')
xlabel('Dirección del viento [º]')
ylabel('Presión [Pa]')
grid on

%% error WS y WD

error_WS = zeros(1,length(datosTotal.ICAO));
error_WD = zeros(1,length(datosTotal.ICAO));

for i=1:length(datosTotal.ICAO)
    error_WS(i)=abs(datosTotal.WindSpeed(i) - datosTotal.WS_era5(i));
    error_WD(i)=abs(datosTotal.WindDirec_(i) - datosTotal.WD_era5(i));
end

%% WS y WD ERA5 por altura
figure
subplot(1,2,1)
boxplot(abs(datosTotal.WindSpeed - datosTotal.WS_era5),Id2,'orientation','horizontal','Labels',0:14)
title('Velocidad del viento ERA5 por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

subplot(1,2,2)
boxplot(abs(datosTotal.WindDirec_ - datosTotal.WD_era5),Id2,'orientation','horizontal','Labels',0:14)
title('Dirección del viento ERA5 por nivel Presión')
xlabel('Dirección del viento [º]')
ylabel('Presión [Pa]')
grid on

%% Media de datos por nivel de presión
Id4 = logical(Id3);
med_T_presion = zeros(1,27);
med_u_presion = zeros(1,27);
med_v_presion = zeros(1,27);
for i=5:25
    med_T_presion(i)=mean(datosTotal.Temperature(Id4(:,i)));
    med_u_presion(i)=mean(datosTotal.u(Id4(:,i)));
    med_v_presion(i)=mean(datosTotal.v(Id4(:,i)));
end
med_T_presion([1,2,3,4,26,27]) = [];
med_u_presion([1,2,3,4,26,27]) = [];
med_v_presion([1,2,3,4,26,27]) = [];

% Errores por nivel de presión
med_errT_presion = zeros(1,27);
med_erru_presion = zeros(1,27);
med_errv_presion = zeros(1,27);
for i=5:25
    med_errT_presion(i)=mean(datosTotal.Temperature(Id4(:,i)))-mean(datosTotal.Temp_era5(Id4(:,i)));
    med_erru_presion(i)=mean(datosTotal.u(Id4(:,i)))-mean(datosTotal.u_era5(Id4(:,i)));
    med_errv_presion(i)=mean(datosTotal.v(Id4(:,i)))-mean(datosTotal.v_era5(Id4(:,i)));
end

med_errT_presion([1,2,3,4,26,27]) = [];
med_erru_presion([1,2,3,4,26,27]) = [];
med_errv_presion([1,2,3,4,26,27]) = [];


%Media error por Nivel de presión:
errorMedio_T_presion = mean(abs(med_errT_presion))
errorMedio_u_presion =mean(abs(med_erru_presion))
errorMedio_v_presion =mean(abs(med_errv_presion))

%% Media de error de datos por presión
figure
subplot(1,3,1)
boxplot((datosTotal.Temperature-datosTotal.Temp_era5),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Error medio de Temperatura por nivel de Presión')
xlabel('Grados [º]')
ylabel('Presión [Pa]')
grid on

subplot(1,3,2)
boxplot((datosTotal.u -datosTotal.u_era5),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Error medio de u por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

subplot(1,3,3)
boxplot((datosTotal.v -datosTotal.v_era5),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Error medio de v por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

%%
figure
subplot(2,3,1)
boxplot((datosTotal.Temperature),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Temperatura derivada por nivel de Presión')
xlabel('Grados [º]')
ylabel('Presión [Pa]')
grid on

subplot(2,3,2)
boxplot((datosTotal.u),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('u derivada por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

subplot(2,3,3)
boxplot((datosTotal.v),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('v derivada por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

subplot(2,3,4)
boxplot((datosTotal.Temp_era5),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('Temperatura ERA5 por nivel de Presión')
xlabel('Grados [º]')
ylabel('Presión [Pa]')
grid on

subplot(2,3,5)
boxplot((datosTotal.u_era5),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('u ERA5 por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on

subplot(2,3,6)
boxplot((datosTotal.v_era5),Id3,'orientation','horizontal','Labels',fliplr({'20000' '22500' '25000' '30000' '35000' '40000' '45000' '50000' '55000' '60000' '65000' '70000' '75000' '77500' '80000' '82500' '85000' '87500' '90000' '92500' '95000'}))
title('v ERA5 por nivel Presión')
xlabel('Velocidad del viento [m/s]')
ylabel('Presión [Pa]')
grid on


