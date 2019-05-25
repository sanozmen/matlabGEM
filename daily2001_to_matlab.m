%a photosynthesis model calibrated by particularly
%regarding the temperature dependence of Asat and phi(quantum yield)


clear all
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