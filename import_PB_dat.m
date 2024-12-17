%Script to import data for Pleasant Bay
% import_PB_dat.m 
% Inputs:
%   1) Spreadsheets with in situ data along with folder information
%   %
% Outputs:
%   1) In situ data table in mat file
% 
% Last revised: 20 June 2024

%% ------------ Begin Code -------------
clc
clearvars

folderpath='C:\Users\slohrenz\Documents\ArcGIS\Projects\PleasantBay\InSitu_Data\';

%% Import data from spreadsheet (2010 - 2019)

opts1=detectImportOptions([folderpath,'PleasantBayWQdbf_2000to2019_import.xlsx']);
opts1.DataRange = 'A2';
opts1.VariableNamesRange = 'A1';
opts1 = setvartype(opts1,{'Date','SampleTime'},{'datetime','double'});
pb_dat=readtable([folderpath,'PleasantBayWQdbf_2000to2019_import.xlsx'],...
    opts1,"Sheet",'PleasantBayWQdbf_2000to2019');

pb_sta_indx=find(contains(pb_dat.Program,'PBA'));
pb_sta=pb_dat.StationID(pb_sta_indx);
pb_date=pb_dat.Date(pb_sta_indx);
pb_year=pb_dat.Year(pb_sta_indx);
pb_time=pb_dat.SampleTime(pb_sta_indx);
pb_totdep=pb_dat.Total_depth_m_(pb_sta_indx);
pb_dep=pb_dat.Depth(pb_sta_indx);
pb_chl=pb_dat.Chla_ug_L_(pb_sta_indx);
pb_tot_pig=pb_dat.Total_Pigments_ug_L_(pb_sta_indx);
pb_secchi=pb_dat.Secchi_avg_m_(pb_sta_indx);
pb_do=pb_dat.FieldDO_conc_mg_L_(pb_sta_indx);

%% Import data from spreadsheet (2020)

opts2=detectImportOptions([folderpath,'CHATHAM 2020_Client Final_edited_for_import.xlsx']);
opts2.DataRange = 'A18';
opts2.VariableNamingRule = 'preserve';
opts2.VariableNamesRange = 'A17';
pb2_dat=readtable([folderpath,'CHATHAM 2020_Client Final_edited_for_import.xlsx'],...
    opts2,"Sheet",'CHATHAM 2020 Client Final');

pb2_sta_indx=find(contains(pb2_dat.SampleID,'PBA'));
pb2_sta_num=pb2_dat.("Station#")(pb2_sta_indx);
for ista = 1:length(pb2_sta_indx)
    pb2_sta{ista}=['PBA-',strrep(num2str(pb2_sta_num(ista),'%u'),' ','')];
end
pb2_sta=pb2_sta';
pb2_date=pb2_dat.Date(pb2_sta_indx);
pb2_time=pb2_dat.SAMPLETIME(pb2_sta_indx);
pb2_totdep=pb2_dat.TOTDEP(pb2_sta_indx);
pb2_dep=pb2_dat.Depth(pb2_sta_indx);
pb2_chl=pb2_dat.CHLMICROGPERL(pb2_sta_indx);
pb2_tot_pig=pb2_dat.TOTPIGMICROGPERL(pb2_sta_indx);
pb2_secchi=pb2_dat.SECCHIAVERAGE(pb2_sta_indx);
pb2_do=pb2_dat.DOMGPERL(pb2_sta_indx);

%% Import data from spreadsheet (2021)

opts3=detectImportOptions([folderpath,'CHATHAM 2021 Client Final_edited_for_import.xlsx']);
opts3.DataRange = 'A18';
opts3.VariableNamingRule = 'preserve';
opts3.VariableNamesRange = 'A17';
pb3_dat=readtable([folderpath,'CHATHAM 2021 Client Final_edited_for_import.xlsx'],opts2,"Sheet",'CLIENT FINAL');

pb3_sta_indx=find(contains(pb3_dat.SampleID,'PBA'));
pb3_sta_num=pb3_dat.("Station#")(pb3_sta_indx);
for ista = 1:length(pb3_sta_indx)
    pb3_sta{ista}=['PBA-',strrep(num2str(pb3_sta_num(ista),'%u'),' ','')];
end
pb3_sta=pb3_sta';
pb3_date=pb3_dat.Date(pb3_sta_indx);
pb3_time=pb3_dat.SAMPLETIME(pb3_sta_indx);
pb3_totdep=pb3_dat.TOTDEP(pb3_sta_indx);
pb3_dep=pb3_dat.Depth(pb3_sta_indx);
pb3_chl=pb3_dat.CHLMICROGPERL(pb3_sta_indx);
pb3_tot_pig=pb3_dat.TOTPIGMICROGPERL(pb3_sta_indx);
pb3_secchi=pb3_dat.SECCHIAVERAGE(pb3_sta_indx);
pb3_do=pb3_dat.DOMGPERL(pb3_sta_indx);

%% Merge data
pb_sta = [pb_sta;pb2_sta;pb3_sta];
pb_date = [pb_date;pb2_date;pb3_date];
pb_time = [pb_time;pb2_time;pb3_time];
pb_date_time = pb_date + days(pb_time) + hours(4);  % Adjust to UTC
pb_totdep = [pb_totdep;pb2_totdep;pb3_totdep];
pb_dep = [pb_dep;pb2_dep;pb3_dep];
pb_chl = [pb_chl;pb2_chl;pb3_chl];
pb_tot_pig = [pb_tot_pig;pb2_tot_pig;pb3_tot_pig];
pb_secchi = [pb_secchi;pb2_secchi;pb3_secchi];
pb_do = [pb_do;pb2_do;pb3_do];

%% Save data

save([folderpath,'all_pb_dat_2010_2021.mat'],'pb_sta','pb_date','pb_time','pb_date_time','pb_totdep','pb_dep',...
    'pb_chl','pb_tot_pig','pb_secchi','pb_do');

disp('Completed');


