clear
close all
Run=1
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


%Data preparation
% Negative GPP is converted to NaN
GPP(GPP<0)=NaN;
PAR(PAR<0)=NaN;
I=PAR(7000:9000); %subset for growing season
I=I(I<200 & I>100);
GPP_s=GPP(7000:9000); % subset for growing season
GPP_s=GPP_s(I<200 & I>100);
%parameters
%irradiance

%temperature
T=AirT(7000:9000);
T=T(I<200 & I>100);
%constant ambient co2 concentration µmol mol^-1
Ca=CO2c(7000:9000);
Ca=Ca(I<200 & I>100);
%beta µmol mol^-1
beta=100;
%phi15 notional max value at 15 centigrat µmol mol^-1
phi15=0.098;
%template for f_tphi
f_tphi=zeros(size(T));
%constant coefficient of carbon dependence of Asat
Kca_Asat=200;
%Temperature dependence of Asat
T_0=-5;
Tref=20;
Tmax_Asat350=28;
Tmax_Asat200=32;

%1 PHI /////////////////
%phi which is the quantum yield
%Carbon dependence of phi
f_caphi=1-beta./Ca;

%temperature dependence of phi
f_tphi=1-(0.015.*(T-15).*(350./Ca));
phi=phi15.*f_caphi.*f_tphi;
phi(T<=15)=phi15.*f_caphi(T<=15);
  
%2 ASAT ////////////////

%equation describing carbon dependence of Asat
M=Kca_Asat./Ca;
f_caAsat=1./(1+M);

Tmax_Asat=Tmax_Asat350+(Tmax_Asat200-Tmax_Asat350).*((Ca-350)./(200-350));

T_0_Asat=(3.*Tmax_Asat-T_0)./(2);

%template for f_TAsat
%calculate value of f_A_sat
f_TAsat=(((T-T_0).^2).*(T_0_Asat-T))./(((Tref-T_0).^2).*(T_0_Asat-Tref));
f_TAsat(T<T_0)=0;


%3 FINAL  ///////////////////

%notional maximum value when other factors are unity at 20 oC (µmolCO2 m-2s-1)
A_sat20=60;
%convexity factor teta
teta=0.5;

        
%equation used to describe the temperature and CO2 tependence of A_sat
A_sat=A_sat20.*f_caAsat.*f_TAsat;


 A_stem=(phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - 4.*teta.*phi.*I.*A_sat))./(2.*teta);


sq_diff=sum((GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem)))).^2); % Calculate the fit of the model to observed data
old_diff=sq_diff;

%Monte Carlo analysis to find optimal values of teta and A_sata20 using the
%lowest square difference as measurement of goodnes of fit 

    for teta=0.1:0.01:1;
                         
                 A_stem=(phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - (4.*teta.*phi.*I.*A_sat)))./(2.*teta);
                 sq_diff=sum((GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem)))).^2);
            if sq_diff<old_diff;
            old_diff=sq_diff;
            optimal_teta=teta;
            optimal_A_sat20=A_sat20;
          
          
    end
end


optimal_teta
optimal_A_sat20
%Model using optimal values

A_sat=optimal_A_sat20.*f_caAsat.*f_TAsat;
A_stem=(phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - 4.*optimal_teta.*phi.*I.*A_sat))./(2.*optimal_teta);


%Statistics of goodness of fit
n=numel(GPP_s(isfinite(GPP_s)& isfinite(A_stem))); %number of elements used in the analysis
O=GPP_s(isfinite(GPP_s)& isfinite(A_stem)); %observations which are not NaN in both vectors
P=A_stem(isfinite(GPP_s)& isfinite(A_stem));%model which are not NaN in both vectors
%correlation coefficient
R2=corrcoef(GPP_s(isfinite(GPP_s)& isfinite(A_stem)), A_stem(isfinite(GPP_s)& isfinite(A_stem)));
RMSE=sqrt(mean((GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem)))).^2));
RMSE_percent=RMSE.*100./mean(GPP_s(isfinite(GPP_s)));
mean_diff=mean(GPP_s(isfinite(GPP_s)& isfinite(A_stem))-(A_stem(isfinite(GPP_s)& isfinite(A_stem))));

T_test=(mean_diff.*sqrt(n))./sqrt(sum((O - P - mean_diff).^2)./(n-1));

%A table with all the statistics 
%R2  RMSE   RMSE%   M   t_test     optimal_teta   Optimal_Asat20 (µmmol) 
Results(Run,1)=R2(1,2);
Results(Run,2)=RMSE;
Results(Run,3)=RMSE_percent;
Results(Run,4)=mean_diff;
Results(Run,5)=T_test
Results(Run,6)=optimal_teta;
Results(Run,7)=optimal_A_sat20;


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



%run 99 data
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

%Data preparation
% Negative GPP is converted to NaN
GPP(GPP<0)=NaN;
PAR(PAR<0)=NaN;
I=PAR(7000:9000); %subset for growing season
I=I(I<200 & I>100);

GPP_s=GPP(7000:9000); % subset for growing season
GPP_s=GPP_s(I<200 & I>100);
%parameters
%irradiance
%temperature
T=AirT(7000:9000);
T=T(I<200 & I>100);
%constant ambient co2 concentration µmol mol^-1
Ca=CO2c(7000:9000);
Ca=Ca(I<200 & I>100);

%1 PHI /////////////////%run 99 data
%run 99 data
%Carbon dependence of phi
f_caphi=1-beta./Ca;

%temperature dependence of phi
f_tphi=1-(0.015.*(T-15).*(350./Ca));
phi=phi15.*f_caphi.*f_tphi;
phi1=phi15.*f_caphi.*1;
phi(T<=15)=phi1(T<=15);

%2 ASAT ////////////////%run 99 data

%constant coefficient of carbon dependence of Asat
Kca_Asat=200;
%equation describing carbon dependence of Asat
M=Kca_Asat./Ca;
f_caAsat=1./(1+M);

T_0_Asat=(3.*Tmax_Asat-T_0)./(2);


  %template for f_TAsat
%calculate value of f_A_sat
f_TAsat=(((T-T_0).^2).*(T_0_Asat-T))./(((Tref-T_0).^2).*(T_0_Asat-Tref));
f_TAsat(T<T_0)=0;


% Run model with 99 data
A_sat=optimal_A_sat20.*f_caAsat.*f_TAsat;
A_stem=(phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - 4.*optimal_teta.*phi.*I.*A_sat))./(2.*optimal_teta);


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
Results(Run,7)=optimal_A_sat20;


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


%run 2001 data
%transfer_csv_data_to_matlab
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
% Negative GPP is converted to NaN
GPP(GPP<0)=NaN;
PAR(PAR<0)=NaN;
I=PAR(7000:9000); %subset for growing season
I=I(I<200 & I>100);

GPP_s=GPP(7000:9000); % subset for growing season
GPP_s=GPP_s(I<200 & I>100);
%parameters
%irradiance
I=PAR(7000:9000); %subset for growing season
I=I(I<200 & I>100);
%temperature
T=AirT(7000:9000);
T=T(I<200 & I>100);
%constant ambient co2 concentration µmol mol^-1
Ca=CO2c(7000:9000);
Ca=Ca(I<200 & I>100);

%1 PHI /////////////////%run 01 data
%run 99 data
%Carbon dependence of phi
f_caphi=1-beta./Ca;

%temperature dependence of phi
f_tphi=1-(0.015.*(T-15).*(350./Ca));
phi=phi15.*f_caphi.*f_tphi;
phi(T<=15)=phi15.*f_caphi(T<=15);
  

%2 ASAT ////////////////%run 01 data

%equation describing carbon dependence of Asat
M=Kca_Asat./Ca;
f_caAsat=1./(1+M);

T_0_Asat=(3.*Tmax_Asat-T_0)./(2);

%template for f_TAsat
%calculate value of f_A_sat
f_TAsat=(((T-T_0).^2).*(T_0_Asat-T))./(((Tref-T_0).^2).*(T_0_Asat-Tref));
f_TAsat(T<T_0)=0;


% Run model with 01 data
A_sat=optimal_A_sat20.*f_caAsat.*f_TAsat;
A_stem=(phi.*I + A_sat - sqrt((phi.*I + A_sat).^2 - 4.*optimal_teta.*phi.*I.*A_sat))./(2.*optimal_teta);




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
Results(Run,7)=optimal_A_sat20;


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
