close all
clear all
%a photosynthesis model using non rectangular calibrated by particularly
%regarding the temperature dependence of Asat and phi(quantum yield)

    %////////////////////////////CALIBRATION ////////////////////////////////////
   %//////////////////////////////////////////////////////////////////////////
  %/////////////////////////////////////////////////////////////////////////
Run=1;% first run of the model with 99 data to calibrate and obtain optimal values for 
%A_sat20 and Teta
%transfer_csv_data_to_matlab
%From 30min data we use only PAR 
fid=fopen('C:\TempData\em\Norunda\30min_data\Norunda2001.csv','r');%here you have to change the folder and file name
headerline=fgetl(fid);%read the first line 
%replace komma by space
headerline=strrep(headerline,',',' ');
%read all data starting from the line number given as the second argument to csvread below:
%NOTE: 0-based, e.g. 2=third line. Note: check file first to find out line number for first time step
all_data = csvread('C:\TempData\em\Norunda\30min_data\Norunda2001.csv',1);
for i=1:size(all_data,2)
    [varname,headerline]=strtok(headerline); %extract variable names
    
    eval([varname '=all_data(:,',num2str(i),');']);% assign the data to the varname
end
%Data preparation
%-9999 data to NaN
PAR(PAR==-9999)=NaN;

PAR=PAR*1800;%PAR is recalculated to �mol/m2x30min
B=reshape(PAR,48,365); %each 48 (30 min) values make one day of the year
Daily_PAR=(sum(B)/1000000)'; %from �mol/m2d -> mol/m2d



%transfer_csv_data_to_matlab
fid=fopen('C:\TempData\em\Norunda\Daily_data\Norunda2001day.csv','r');%here you have to change the folder and file name
headerline=fgetl(fid);%read the first line 
%replace komma by space
headerline=strrep(headerline,',',' ');
%read all data starting from the line number given as the second argument to csvread below:
%NOTE: 0-based, e.g. 2=third line. Note: check file first to find out line number for first time step
all_data = csvread('C:\TempData\em\Norunda\Daily_data\Norunda2001day.csv',1);
for i=1:size(all_data,2)
    [varname,headerline]=strtok(headerline); %extract variable names
    
    eval([varname '=all_data(:,',num2str(i),');']);% assign the data to the varname
end


%Data preparation
% Negative GPP is converted to NaN
GPP(GPP<0)=NaN;
GPP_s=GPP(90:280); % subset for growing season
%parameters
%irradiance
I=Daily_PAR(90:280); %subset for growing season
%temperature
T=AirT(90:280);
%constant ambient co2 concentration �mol mol^-1
Ca=370;
%beta �mol mol^-1
beta=100;
%phi15 notional max value at 15 centigrat �mol mol^-1
phi15=0.098;
%template for f_tphi
f_tphi=zeros(size(T));
%constant coefficient of carbon dependence of Asat
Kca_Asat=700;
%Temperature dependence of Asat
T_0=0;
Tref=22;
Tmax_Asat350=28;
Tmax_Asat700=32;

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

%equation describing carbon dependence of Asat
f_caAsat=1/(1+Kca_Asat/Ca);

Tmax_Asat=Tmax_Asat350+(Tmax_Asat700-Tmax_Asat350).*((Ca-350)./(700-350));

T_0_Asat=(3.*Tmax_Asat-T_0)./(2);

%template for f_TAsat
f_TAsat=zeros(size(T));
%calculate value of f_A_sat
for m=1:length(T)
    k=T(m);
if k<T_0_Asat && k>T_0
    f_TAsat(m,:)=(((k-T_0).^2).*(T_0_Asat-k))./(((Tref-T_0).^2).*(T_0_Asat-Tref));
else
   f_TAsat(m,:)=0;
end
end


%3 FINAL  ///////////////////

%notional maximum value when other factors are unity at 20 oC (�molCO2 m-2s-1)
A_sat20inmmol=15;
A_sat20=A_sat20inmmol.*86400./1000000; % (molC02 m-2d-1)
%convexity factor teta
teta=0.01;

        
%equation used to describe the temperature and CO2 tependence of A_sat
A_sat=A_sat20.*f_caAsat.*f_TAsat;


 A_stem=(phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - 4.*teta.*phi.*I.*A_sat))./2.*teta.*12;


sq_diff=sum((GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem)))).^2); % Calculate the fit of the model to observed data
old_diff=sq_diff;

%Monte Carlo analysis to find optimal values of teta and A_sata20 using the
%lowest square difference as measurement of goodnes of fit 
for teta=0.01:0.01:1;
    for A_sat20=15.*86400/1000000:0.001:200.*86400/1000000;
                    if Daily_PAR==NaN
                 A_stem=NaN
              else
                 A_sat=A_sat20.*f_caAsat.*f_TAsat;
                 A_stem=((phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - 4.*teta.*phi.*I.*A_sat))./2.*teta).*12;
                 sq_diff=sum((GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem)))).^2);
            if sq_diff<old_diff;
            old_diff=sq_diff;
            optimal_teta=teta;
            optimal_A_sat20=A_sat20;
                end
                 end
    end
end


optimal_teta
optimal_A_sat20mmol=optimal_A_sat20*1000000/86400

%Model using optimal values
A_sat=optimal_A_sat20.*f_caAsat.*f_TAsat;
A_stem=((phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - 4.*optimal_teta.*phi.*I.*A_sat))./2.*optimal_teta).*12;


%Statistics of goodness of fit
n=numel(GPP_s(isfinite(GPP_s)& isfinite(A_stem))); %number of elements used in the analysis
O=GPP_s(isfinite(GPP_s)& isfinite(A_stem)); %observations which are not NaN in both vectors
P=A_stem(isfinite(GPP_s)& isfinite(A_stem));%model which are not NaN in both vectors
%correlation coefficient
R2=corrcoef(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)));
RMSE=sqrt(mean((GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem)))).^2));
RMSE_percent=RMSE.*100./mean(GPP_s(isfinite(GPP_s)));
mean_diff=mean(GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem))));

T_test=(mean_diff.*sqrt(n))./sqrt(sum((O - P - mean_diff).^2)./(n-1))

%A table with all the statistics 
%R2  RMSE   RMSE%   M   t_test     optimal_teta   Optimal_Asat20 (�mmol) 
Results(Run,1)=R2(1,2);
Results(Run,2)=RMSE;
Results(Run,3)=RMSE_percent;
Results(Run,4)=mean_diff;
Results(Run,5)=T_test
Results(Run,6)=optimal_teta;
Results(Run,7)=optimal_A_sat20mmol;


%Figures 
figure(1) %Plot the model and observation against time to see behavior
plot(A_stem, 'k--')
hold on
plot(GPP_s, 'b')
hold off

figure(2)
plot(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)), 'o')
hold on
plot(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), polyval(polyfit(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)),1),GPP_s(isfinite(GPP_s)& isfinite(A_stem))))
plot([0:14], [0:14], 'k--')
hold off

figure(3)
plot(T, A_stem, 'x')
hold on
plot(T,GPP_s, 'rx')
hold off









%////////// 	EVALUATION   ////////////////////////////////////
%///////////////////////////////////////////////////////
%/////////////////////////////////////////////////
Run=2;% run of the model with 97 data and evaluate using 97 
      %   data using parameters obtained in previuos CALIBRATION run

%transfer_csv_data_to_matlab
fid=fopen('C:\TempData\em\Norunda\30min_data\Norunda1997.csv','r');%here you have to change the folder and file name
headerline=fgetl(fid);%read the first line 
%replace komma by space
headerline=strrep(headerline,',',' ');
%read all data starting from the line number given as the second argument to csvread below:
%NOTE: 0-based, e.g. 2=third line. Note: check file first to find out line number for first time step
all_data = csvread('C:\TempData\em\Norunda\30min_data\Norunda1997.csv',1);
for i=1:size(all_data,2)
    [varname,headerline]=strtok(headerline); %extract variable names
    
    eval([varname '=all_data(:,',num2str(i),');']);% assign the data to the varname
end
%run 99 data
%-9999 data to NaN
PAR(PAR==-9999)=NaN;
PAR=PAR*1800;
B=reshape(PAR,48,365);
Daily_PAR=(sum(B)/1000000)';

%run 99 data
%transfer_csv_data_to_matlab
fid=fopen('C:\TempData\em\Norunda\Daily_data\Norunda1997day.csv','r');%here you have to change the folder and file name
headerline=fgetl(fid);%read the first line 
%replace komma by space
headerline=strrep(headerline,',',' ');
%read all data starting from the line number given as the second argument to csvread below:
%NOTE: 0-based, e.g. 2=third line. Note: check file first to find out line number for first time step
all_data = csvread('C:\TempData\em\Norunda\Daily_data\Norunda1997day.csv',1);
for i=1:size(all_data,2)
    [varname,headerline]=strtok(headerline); %extract variable names
    
    eval([varname '=all_data(:,',num2str(i),');']);% assign the data to the varname
end

%run 99 data
% Negative GPP is converted to NaN
GPP(GPP<0)=NaN;
GPP_s=GPP(90:280); % subset for growing season
%parameters
%irradiance
I=Daily_PAR(90:280); %subset for growing season
%temperature
T=AirT(90:280);
%constant ambient co2 concentration �mol mol^-1
Ca=370;
%beta �mol mol^-1
beta=100;
%phi15 notional max value at 15 centigrat �mol mol^-1
phi15=0.098;
%template for f_tphi
f_tphi=zeros(size(T));

%1 PHI /////////////////%run 99 data
%run 99 data
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

%2 ASAT ////////////////%run 99 data

%constant coefficient of carbon dependence of Asat
Kca_Asat=700;
%equation describing carbon dependence of Asat
f_caAsat=1/(1+Kca_Asat/Ca);

T_0_Asat=(3.*Tmax_Asat-T_0)./(2);

%template for f_TAsat
f_TAsat=zeros(size(T));

for m=1:length(T)
    k=T(m);
if k<T_0_Asat && k>T_0
    f_TAsat(m,:)=(((k-T_0).^2).*(T_0_Asat-k))./(((Tref-T_0).^2).*(T_0_Asat-Tref));
else
   f_TAsat(m,:)=0;
end
end

% Run model with 99 data
A_sat=optimal_A_sat20.*f_caAsat.*f_TAsat;
A_stem=((phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - 4.*optimal_teta.*phi.*I.*A_sat))./2.*optimal_teta).*12;


%run 99 data
% Evaluation of the model
n=numel(GPP_s(isfinite(GPP_s)& isfinite(A_stem))); %number of elements used in the analysis
O=GPP_s(isfinite(GPP_s)& isfinite(A_stem)); %observations which are not NaN in both vectors
P=A_stem(isfinite(GPP_s)& isfinite(A_stem));%model which are not NaN in both vectors
%correlation coefficient
R2=corrcoef(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)));
RMSE=sqrt(mean((GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem)))).^2));
RMSE_percent=RMSE.*100./mean(GPP_s(isfinite(GPP_s)));
mean_diff=mean(GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem))));

T_test=(mean_diff.*sqrt(n))./sqrt(sum((O - P - mean_diff).^2)./(n-1));

Results(Run,1)=R2(1,2);
Results(Run,2)=RMSE;
Results(Run,3)=RMSE_percent;
Results(Run,4)=mean_diff;
Results(Run,5)=T_test
Results(Run,6)=optimal_teta;
Results(Run,7)=optimal_A_sat20mmol;


%Figures 
%Run 99
figure(Run*10+1) %Plot the model and observation against time to see behavior
plot(A_stem, 'k--')
hold on
plot(GPP_s, 'b')
hold off

figure(Run*10+2)
plot(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)), 'o')
hold on
plot(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), polyval(polyfit(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)),1),GPP_s(isfinite(GPP_s)& isfinite(A_stem))))
plot([0:14], [0:14], 'k--')
hold off

figure(Run*10+3)
plot(T, A_stem, 'x')
hold on
plot(T,GPP_s, 'rx')
hold off





Run=3;% run of the model with 2001 data and evaluate using 01 
      %   data using parameters obtained in previuos CALIBRATION run

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
%run 2001 data
%-9999 data to NaN
PAR(PAR==-9999)=NaN;
PAR=PAR*1800;
B=reshape(PAR,48,365);
Daily_PAR=(sum(B)/1000000)';

%run 2001 data
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


%run 2001 data
% Negative GPP is converted to NaN
GPP(GPP<0)=NaN;
GPP_s=GPP(90:280); % subset for growing season
%parameters
%irradiance
I=Daily_PAR(90:280); %subset for growing season
%temperature
T=AirT(90:280);
%constant ambient co2 concentration �mol mol^-1
%template for f_tphi
f_tphi=zeros(size(T));


%1 PHI /////////////////%run 01 data
%run 99 data
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

%2 ASAT ////////////////%run 01 data

%equation describing carbon dependence of Asat
f_caAsat=1/(1+Kca_Asat/Ca);

T_0_Asat=(3.*Tmax_Asat-T_0)./(2);

%template for f_TAsat
f_TAsat=zeros(size(T));

for m=1:length(T)
    k=T(m);
if k<T_0_Asat && k>T_0
    f_TAsat(m,:)=(((k-T_0).^2).*(T_0_Asat-k))./(((Tref-T_0).^2).*(T_0_Asat-Tref));
else
   f_TAsat(m,:)=0;
end
end

% Run model with 01 data
A_sat=optimal_A_sat20.*f_caAsat.*f_TAsat;
A_stem=((phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - 4.*optimal_teta.*phi.*I.*A_sat))./2.*optimal_teta).*12;




%run 01 data
% Evaluation of the model
n=numel(GPP_s(isfinite(GPP_s)& isfinite(A_stem))); %number of elements used in the analysis
O=GPP_s(isfinite(GPP_s)& isfinite(A_stem)); %observations which are not NaN in both vectors
P=A_stem(isfinite(GPP_s)& isfinite(A_stem));%model which are not NaN in both vectors
%correlation coefficient
R2=corrcoef(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)));
RMSE=sqrt(mean((GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem)))).^2));
RMSE_percent=RMSE.*100./mean(GPP_s(isfinite(GPP_s)));
mean_diff=mean(GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem))));

T_test=(mean_diff.*sqrt(n))./sqrt(sum((O - P - mean_diff).^2)./(n-1));
%Table with statistics
Results(Run,1)=R2(1,2);
Results(Run,2)=RMSE;
Results(Run,3)=RMSE_percent;
Results(Run,4)=mean_diff;
Results(Run,5)=T_test
Results(Run,6)=optimal_teta;
Results(Run,7)=optimal_A_sat20mmol;


%Figures 
%Run 01
figure(Run*10+1) %Plot the model and observation against time to see behavior
plot(A_stem, 'k--')
hold on
plot(GPP_s, 'b')
hold off

figure(Run*10+2)
plot(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)), 'o')
hold on
plot(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), polyval(polyfit(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)),1),GPP_s(isfinite(GPP_s)& isfinite(A_stem))))
plot([0:14], [0:14], 'k--')
hold off

figure(Run*10+3)
plot(T, A_stem, 'x')
hold on
plot(T,GPP_s, 'rx')
hold off
