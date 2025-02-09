% read_matchup_PB.m
% Syntax: read_matchup_PB
%
% Script reads satellite matchup output and compiles into a single .mat
% file
%
% Inputs:
%    1) Directory location for matchup files
%    2) Directory location for output files
%
% Outputs:
%    1) .mat file with compiled satellite matchup data
%   
% Other m-files required: None 
%
% MAT-files required: 
%    1) None
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 26 Sep 2024

%% ------------- BEGIN CODE --------------%

clc
clearvars

% C2RCC
InDir='~\Satellite_matchups\Matchup files\';
OutDir='~\Satellite_matchups\';

% C2X
% InDir = '~\Satellite_matchups\Matchup files\C2X\';
% OutDir='~\Satellite_matchups\Matchup files\C2X\';

% L2gen
% InDir = '~\Satellite_matchups\Matchup files\l2gen\';
% OutDir='~\Satellite_matchups\Matchup files\l2gen\';

% Get list of files

Dlist=dir([InDir,'20*']);
% Dlist=dir([InDir,'\*NOreproject.csv']);  % L2gen
Flist=cell(length(Dlist),1);
fcount = 1;

%% Loop through year directories for all_satup file names
for idir = 1:length({Dlist.name})
    fname=dir([InDir,Dlist(idir).name,'\*.csv']);
    % fname=dir([InDir,Dlist(idir).name]);  % L2gen
    for ik=1:length(fname)
        Flist{fcount}=[fname(ik).folder,'\',fname(ik).name];
        fcount = fcount+1;
    end
end
    
fnum = length(Flist);
ftabl=cell(fnum,1);

%% Define file import options
opts = delimitedTextImportOptions('Datalines',8,'VariableNamesLine',7,'NumVariables',121,...
    'VariableNamingRule','modify');
%opts = setvartype(opts,{'ProdID','CoordID','Name','Latitude','Longitude','Date_yyyy_MM_dd_','Time_HH_mm_ss_'},...
%    {'char','char','double','double','datetime','duration'});
%opts = spreadsheetImportOptions('DataRange','A8','VariableNamesRange','A7','NumVariables',97);
opts = setvartype(opts,{'Var1','Var2','Var3','Var4','Var5','Var6','Var7','Var8','Var9'},...
    {'char','char','char','double','double','double','double','datetime','duration'});
    
%% Read variables for each matchup file for a given year
for ifile=1:fnum
    disp(['Reading ',Flist{ifile}]);
    w=warning('off','MATLAB:table:ModifiedAndSavedVarnames'); 
    ftabl{ifile}=readtable(Flist{ifile},opts);
    warning(w);
    ftabl{ifile}.date_time=ftabl{ifile}.Date_yyyy_MM_dd_ + ftabl{ifile}.Time_HH_mm_ss_;
    flds=fieldnames(ftabl{ifile});
    for ifld=12:53 %length(flds)
        % C2RCC and C2X
        % ftabl{ifile}.(flds{ifld})=str2num(char(ftabl{ifile}.(flds{ifld})));     
        
        % L2gen
        data_out = ftabl{ifile}.(flds{ifld});
        [data_out{strcmp(data_out,'')}] = deal('0'); % Change null values to '0'
        ftabl{ifile}.(flds{ifld})=str2num(char(data_out));
    end
end
    
%% Apportion output by date and station in structured variable
stalist=unique(ftabl{1}.CoordID);
nsta=length(stalist);
all_sta={};

all_sta=struct();
for ifile=1:fnum
    flds=fieldnames(ftabl{ifile});
    flds=flds(~strcmp(flds,'Properties'));  %Eliminate this field
    for ista=1:nsta
        staindx=find(strcmp(stalist(ista),ftabl{ifile}.CoordID));  %CoorID is the station number
        for ifld=[1:53,122] % Includes only fields with data and date_time field (122)
            all_sta(ifile,ista).(flds{ifld})=ftabl{ifile}.(flds{ifld})(staindx);
        end
    end
end

save([OutDir,'all_station_matchup.mat'],'all_sta');

%% Restructure file by date and station

sat_fields = fieldnames(all_sta);
[mdt,nsta]=size(all_sta);  % Get number of image (different dates) and number of stations in all_sta


% C2RCC and C2X
% pb_all_sat = struct('sat_sta',{},'sat_lat',{},'sat_lon',{},'sat_date_time',{},'sat_chl',{},'sat_kd489',{},...
%     'sat_rhow_1',{},'sat_rhow_2',{},'sat_rhow_3',{},'sat_rhow_4',{},'sat_rhow_5',{},'sat_tsm',{});

% L2gen
pb_all_sat = struct('sat_sta',{},'sat_lat',{},'sat_lon',{},'sat_date_time',{},'sat_chl',{},'sat_kd489',{},...
    'sat_rhow_1',{},'sat_rhow_2',{},'sat_rhow_3',{},'sat_rhow_4',{});

count = 1;
for idt=1:mdt
    for ista = 1:nsta
        if count == 1
            pb_all_sat(idt).sat_sta=all_sta(idt,ista).CoordID;
            pb_all_sat(idt).sat_lat=all_sta(idt,ista).Latitude;
            pb_all_sat(idt).sat_lon=all_sta(idt,ista).Longitude;
            pb_all_sat(idt).sat_date_time=all_sta(idt,ista).date_time;
            
            
            % C2RCC and C2X
            pb_all_sat(idt).sat_chl=all_sta(idt,ista).conc_chl;
            pb_all_sat(idt).sat_kd489=all_sta(idt,ista).kd489;
            pb_all_sat(idt).sat_kdz90=all_sta(idt,ista).kd_z90max;
            pb_all_sat(idt).sat_rhow_1=all_sta(idt,ista).rhow_1;
            pb_all_sat(idt).sat_rhow_2=all_sta(idt,ista).rhow_2;
            pb_all_sat(idt).sat_rhow_3=all_sta(idt,ista).rhow_3;
            pb_all_sat(idt).sat_rhow_4=all_sta(idt,ista).rhow_4;
            pb_all_sat(idt).sat_rhow_5=all_sta(idt,ista).rhow_5;
            pb_all_sat(idt).sat_tsm=all_sta(idt,ista).conc_tsm;
 
            % L2gen
            % pb_all_sat(idt).sat_chl=all_sta(idt,ista).chlor_a;
            % pb_all_sat(idt).sat_kd489=all_sta(idt,ista).Kd_490;
            % pb_all_sat(idt).sat_rhow_1=all_sta(idt,ista).Rrs_443;
            % pb_all_sat(idt).sat_rhow_2=all_sta(idt,ista).Rrs_482;
            % pb_all_sat(idt).sat_rhow_3=all_sta(idt,ista).Rrs_561;
            % pb_all_sat(idt).sat_rhow_4=all_sta(idt,ista).Rrs_655;
 
        else
            pb_all_sat.sat_sta=vertcat(pb_all_sat.sat_sta,all_sta(idt,ista).CoordID);
            pb_all_sat.sat_lat=vertcat(pb_all_sat.sat_lat,all_sta(idt,ista).Latitude);
            pb_all_sat.sat_lon=vertcat(pb_all_sat.sat_lon,all_sta(idt,ista).Longitude);
            pb_all_sat.sat_date_time=vertcat(pb_all_sat.sat_date_time,all_sta(idt,ista).date_time);
            
            
            % C2RCC and C2X
            pb_all_sat.sat_chl=vertcat(pb_all_sat.sat_chl,all_sta(idt,ista).conc_chl);
            pb_all_sat.sat_kd489=vertcat(pb_all_sat.sat_kd489,all_sta(idt,ista).kd489);
            pb_all_sat.sat_kdz90=vertcat(pb_all_sat.sat_kdz90,all_sta(idt,ista).kd_z90max);
            pb_all_sat.sat_rhow_1=vertcat(pb_all_sat.sat_rhow_1,all_sta(idt,ista).rhow_1);
            pb_all_sat.sat_rhow_2=vertcat(pb_all_sat.sat_rhow_2,all_sta(idt,ista).rhow_2);
            pb_all_sat.sat_rhow_3=vertcat(pb_all_sat.sat_rhow_3,all_sta(idt,ista).rhow_3);
            pb_all_sat.sat_rhow_4=vertcat(pb_all_sat.sat_rhow_4,all_sta(idt,ista).rhow_4);
            pb_all_sat.sat_rhow_5=vertcat(pb_all_sat.sat_rhow_5,all_sta(idt,ista).rhow_5);
            pb_all_sat.sat_tsm=vertcat(pb_all_sat.sat_tsm,all_sta(idt,ista).conc_tsm);
            
            % L2gen
            % pb_all_sat.sat_chl=vertcat(pb_all_sat.sat_chl,all_sta(idt,ista).chlor_a);
            % pb_all_sat.sat_kd489=vertcat(pb_all_sat.sat_kd489,all_sta(idt,ista).Kd_490);
            % pb_all_sat.sat_rhow_1=vertcat(pb_all_sat.sat_rhow_1,all_sta(idt,ista).Rrs_443);
            % pb_all_sat.sat_rhow_2=vertcat(pb_all_sat.sat_rhow_2,all_sta(idt,ista).Rrs_482);
            % pb_all_sat.sat_rhow_3=vertcat(pb_all_sat.sat_rhow_3,all_sta(idt,ista).Rrs_561);
            % pb_all_sat.sat_rhow_4=vertcat(pb_all_sat.sat_rhow_4,all_sta(idt,ista).Rrs_655);
        end
        count = count + 1;
    end
end

save([OutDir,'pb_all_sat.mat'],'pb_all_sat');

%{
mean_chl=zeros(nsta,fnum);
for ifile=1:fnum
    for ista=1:nsta
        staindx=find(strcmp(stalist(ista),ftabl{ifile}.Name));
        mean_chl(ista,ifile)=mean(ftabl{ifile}.conc_chl(staindx));
        mean_kd489(ista,ifile)=mean(ftabl{ifile}.kd489(staindx));
        mean_tsm(ista,ifile)=mean(ftabl{ifile}.conc_tsm(staindx));
    end
end
%}

disp('Completed');
