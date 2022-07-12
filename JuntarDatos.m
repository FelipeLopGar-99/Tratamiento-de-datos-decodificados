clear
clc

% Lectura de ficheros
delimeters = {'tab','\t'};
T08 = readtable('Datos_Decodificados_08.txt','Delimiter',delimeters);
T12 = readtable('Datos_Decodificados_12.txt','Delimiter',delimeters);

% Eliminación de Atípicos
datAtipic08 = [];
datAtipic12 = [];
for i=1:(length(T08.Alt))
    if T08.WindSpeed(i) > 111.1 | T08.Alt(i) > 15000.0 | T08.Roll(i) > abs(5.0) | T08.Mach(i) > 1.0 | T08.Temperature(i) > 320
        datAtipic08 = [datAtipic08 i];
    end
    if mod(i,40000) == 0
        i
    end
end
T08(datAtipic08,:) = [];

for i=1:(length(T12.Alt))
    if T12.WindSpeed(i) > 111.1 | T12.Alt(i) > 15000.0 | abs(T12.Roll(i)) > 5.0 | T12.Mach(i) > 1.0 | T12.Temperature(i) > 320
        datAtipic12 = [datAtipic12 i];
    end
    if mod(i,40000) == 0
        i
    end
end
T12(datAtipic12,:) = [];

% Concatenación de datos
T = [T08;T12];
fprintf('Fin \n');


writetable(T,'DatosCompletos.txt','Delimiter','\t','WriteRowNames',true);
type DatosCompletos.txt;

writetable(T08,'Datos08.txt','Delimiter','\t','WriteRowNames',true);
type Datos08.txt;

writetable(T12,'Datos12.txt','Delimiter','\t','WriteRowNames',true);
type Datos12.txt;


%% Juntar datos BDS44:

format longG
% Lectura de ficheros
delimeters = {'tab','\t'};
T08_44 = readtable('BDS44_08.txt','Delimiter',delimeters,'Format','%s%s%s%s%s%s%s%s');
T12_44 = readtable('BDS44_12.txt','Delimiter',delimeters,'Format','%s%s%s%s%s%s%s%s');

% Eliminación de ceros
ceros08 = [];
ceros12 = [];
for i=1:(length(T08_44.TimeStamp))
    if ismember(T08_44.WindSpeed(i),'NONE') | ismember(T08_44.WindDirec(i),'NONE') | ismember(T08_44.Temperature(i),'NONE')
        ceros08 = [ceros08 i];
    end
end
T08_44(ceros08,:) = [];

for i=1:(length(T12_44.TimeStamp))
    if ismember(T12_44.WindSpeed(i),'NONE') | ismember(T12_44.WindDirec(i),'NONE') | ismember(T12_44.Temperature(i),'NONE')
        ceros12 = [ceros12 i];
    end
end
T12_44(ceros12,:) = [];

T_44 = [T08_44;T12_44];

%Quitamos datos de presión, humedad y turbulencia por que no tenemos ningún valor
T_44 = removevars(T_44,{'Pressure','Humidity', 'Turbulence'});

%Convertimos el tipo de datos al deseado
T_44.TimeStamp = str2double(T_44.TimeStamp);
T_44.WindSpeed = str2double(T_44.WindSpeed);
T_44.WindDirec = str2double(T_44.WindDirec);
T_44.Temperature = 273.15 + str2double(T_44.Temperature);

writetable(T_44,'DatosCompletos_BDS44.txt','Delimiter','\t','WriteRowNames',true);
type DatosCompletos_BDS44.txt;
