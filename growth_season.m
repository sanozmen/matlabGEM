close all
clear all
%a photosynthesis model calibrated by particularly
%regarding the temperature dependence of Asat and phi(quantum yield)



%transfer_csv_data_to_matlab
fid=fopen('C:\TempData\em\Norunda\30min_data\Norunda1999.csv','r');%here you have to change the folder and file name
headerline=fgetl(fid);%read the first line 
%replace komma by space
headerline=strrep(headerline,',',' ');
%read all data starting from the line number given as the second argument to csvread below:
%NOTE: 0-based, e.g. 2=third line. Note: check file first to find out line number for first time step
all_data = csvread('C:\TempData\em\Norunda\30min_data\Norunda1999.csv',1);
for i=1:size(all_data,2)
    [varname,headerline]=strtok(headerline); %extract variable names
    
    eval([varname '=all_data(:,',num2str(i),');']);% assign the data to the varname
end

%-9999 data to NaN
PAR(PAR==-9999)=0;
%/////////////////////////////////////////////////////
B=reshape(PAR,48,365);
Daily_PAR=(sum(B)/48)';
%///////////////////////////////////////////////////

Growing_PAR=Daily_PAR(90:275);

%transfer_csv_data_to_matlab
fid=fopen('C:\TempData\em\Norunda\Daily_data\Norunda1999day.csv','r');%here you have to change the folder and file name
headerline=fgetl(fid);%read the first line 
%replace komma by space
headerline=strrep(headerline,',',' ');
%read all data starting from the line number given as the second argument to csvread below:
%NOTE: 0-based, e.g. 2=third line. Note: check file first to find out line number for first time step
all_data = csvread('C:\TempData\em\Norunda\Daily_data\Norunda1999day.csv',1);
for i=1:size(all_data,2)
    [varname,headerline]=strtok(headerline); %extract variable names
    
    eval([varname '=all_data(:,',num2str(i),');']);% assign the data to the varname
end



%-9999 data to NaN
PAR(PAR==-9999)=0;
% Negative GPP is converted to NaN
GPP(GPP<0)=0;
%parameters
%irradiance
I=PAR(:,1);
%temperature
T=AirT(:,1);
%constant ambient co2 concentration �mol mol^-1
Ca=370;
%beta �mol mol^-1
beta=100;
%phi15 notional max value at 15 centigrat �mol mol^-1
phi15=0.098;
%template for f_tphi
f_tphi=zeros(size(T));


%1 PHI /////////////////

%Carbon dependence of phi
f_caphi=1-beta./Ca;

%temperature dependence of phi
for i=1:length(T)
    t=T(i);
if t<=15
    f_tphi(i,:)=1;
else
    f_tphi(i,:)=1-0.015.*(350./Ca);
end
end

%phi which is the quantum yield
phi=phi15.*f_caphi.*f_tphi;

%2 ASAT ////////////////

%constant coefficient of carbon dependence of Asat
Kca_Asat=700;
%equation describing carbon dependence of Asat
f_caAsat=1/(1+Kca_Asat/Ca);

%Temperature dependence of Asat
T_0=0;
Tref=20;
Tmax_Asat350=28;
Tmax_Asat700=30;
Tmax_Asat=Tmax_Asat350+(Tmax_Asat700-Tmax_Asat350).*((Ca-350)./(700-350));

T_0_Asat=(3.*Tmax_Asat-T_0)./(2);

%template for f_TAsat
f_TAsat=zeros(size(T));

for m=1:length(T)
    k=T(m);
if k<T_0_Asat && k>T_0
    f_TAsat(m,:)=((k-T_0).^2.*(T_0_Asat-k))./((Tref-T_0).^2.*(T_0_Asat-Tref));
else
   f_TAsat(m,:)=0;
end
end


%3 FINAL  ///////////////////

%notional maximum value when are unity at 20 oC (�molCO2 m-2s-1)
A_sat20=15;
%convexity factor teta
teta=0.01;

%equation used to describe the temperature and CO2 tependence of A_sat
A_sat=A_sat20.*f_caAsat.*f_TAsat;

% Non rectangular hyperbola model of gross photosynthesis (�mol CO2 m-2)
A_stem=(phi(90:275).*Daily_PAR(90:275) + A_sat(90:275) - sqrt((phi(90:275).*Daily_PAR(90:275) + A_sat(90:275)).^2 - 4.*teta.*phi(90:275).*Daily_PAR(90:275).*A_sat(90:275)))./2.*teta;
Distance=sum((GPP(90:275)-A_stem).^2); % Calculate the fit of the model to observed data
old_distance=Distance;

for teta=0.01:0.01:1;
    for A_sat20=15:1:45;  
        A_sat=A_sat20.*f_caAsat.*f_TAsat;
        A_stem=(phi(90:275).*Daily_PAR(90:275) + A_sat(90:275) - sqrt((phi(90:275).*Daily_PAR(90:275) + A_sat(90:275)).^2 - 4.*teta.*phi(90:275).*Daily_PAR(90:275).*A_sat(90:275)))./2.*teta;
        Distance=sum((GPP(90:275)-A_stem).^2);
        if Distance<old_distance;
            old_distance=Distance;
            optimal_teta=teta;
            optimal_A_sat20=A_sat20;
        end
    end
end

optimal_teta
optimal_A_sat20


plot(A_stem)
hold on
plot(GPP(90:275), 'rx')
hold off



