% matchup_PB_sat_insitu.m
% Syntax: matchup_PB_sat_insitu
%
% Script loads in_situ data and satellite matchup data files and finds
% colocated data within a designated matchup window
%
% Inputs:
%    1) Directory location for satellite matchup .mat file for each station
%    and satellite date generated with read_matchup_PB.m
%    2) Directory location for in situ .mat file generated with import_PB.m
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
% Last revision: 17 June 2024

%% ------------- BEGIN CODE --------------%

clc
clearvars

%% Load data
 
% Select algorithm 

algthm = 'c2rcc';

switch algthm
    case 'c2rcc'
        % C2RCC
        inputarray='~\InSitu_Data\all_pb_dat_2010_2021.mat';  %For evaluating algorithm on single day L3 files
        matcharray='~\Satellite_matchups\all_station_matchup.mat';
        outdir='~\Satellite_matchups\';
    case 'c2x'
        % C2X
        inputarray='~\InSitu_Data\all_pb_dat_2010_2021.mat';  %For evaluating algorithm on single day L3 files
        matcharray='~\Satellite_matchups\Matchup files\C2X\all_station_matchup.mat';
        outdir='~\Satellite_matchups\Matchup files\C2X\';
    case 'l2gen'
        % L2gen
        inputarray='~\InSitu_Data\all_pb_dat_2010_2021.mat';  %For evaluating algorithm on single day L3 files
        matcharray='~\Satellite_matchups\Matchup files\l2gen\all_station_matchup.mat';
        outdir='~\Satellite_matchups\Matchup files\l2gen\';
end

all_pb_match=[];

% Load in situ data
load(inputarray);

% Select matchup window
hour_window = 36; % 36

% Variable list:  'pb_sta','pb_date','pb_time','pb_date_time','pb_totdep','pb_dep',
% 'pb_chl','pb_tot_pig','pb_secchi','pb_do'

% Load satellite matchup data 
load(matcharray);
sat_fields = fieldnames(all_sta);

%% Compare in situ results by image date 

[mdt,nsta]=size(all_sta);  % Get number of image (different dates) and number of stations in all_sta

% Find matching dates (and times?) and pixel locations
switch algthm
    case 'c2rcc'
        % C2RCC and C2X
        pb_match=struct('insitu_sta',{},'sat_sta',{},'sat_lat',{},'sat_lon',{},'sat_date_time',{},'insitu_date_time',{},...
            'insitu_totdepth',{},'insitu_chl',{},'sat_chl',{},'insitu_secchi',{},'sat_kd489',{},'sat_kdz90',{},...
            'insitu_do',{},'sat_rhow_1',{},'sat_rhow_2',{},'sat_rhow_3',{},'sat_rhow_4',{},'sat_rhow_5',{},'sat_tsm',{});
    case 'c2x'
        % C2RCC and C2X
        pb_match=struct('insitu_sta',{},'sat_sta',{},'sat_lat',{},'sat_lon',{},'sat_date_time',{},'insitu_date_time',{},...
            'insitu_totdepth',{},'insitu_chl',{},'sat_chl',{},'insitu_secchi',{},'sat_kd489',{},'sat_kdz90',{},...
            'insitu_do',{},'sat_rhow_1',{},'sat_rhow_2',{},'sat_rhow_3',{},'sat_rhow_4',{},'sat_rhow_5',{},'sat_tsm',{});
    case 'l2gen'
        % L2gen
        pb_match=struct('insitu_sta',{},'sat_sta',{},'sat_lat',{},'sat_lon',{},'sat_date_time',{},'insitu_date_time',{},...
            'insitu_totdepth',{},'insitu_chl',{},'sat_chl',{},'insitu_secchi',{},'sat_kd489',{},...
            'insitu_do',{},'sat_rhow_1',{},'sat_rhow_2',{},'sat_rhow_3',{},'sat_rhow_4',{});  % C2RCC and C2X ,'sat_rhow_5',{},'sat_tsm',{});
end

day_count=0;
for idt=1:mdt
    dt_indx = find(abs(pb_date_time-all_sta(idt,1).date_time(1)) < hours(hour_window) & ...
        contains(pb_dep,'S') & all(~strcmp(all_sta(idt,1).CoordID,'2')));
    if isempty(dt_indx)
        continue
    end
    day_count = day_count + 1;% Match results by station

    for isamp=1:length(dt_indx)
        pb_match(day_count,isamp).insitu_sta=pb_sta(dt_indx(isamp));
        sta_num=char(pb_sta(dt_indx(isamp)));
        sta_num=(strrep(sta_num(5:end),'A',''));
        sat_sta_indx = [];
        for ista=1:21
            if strcmp(all_sta(idt,ista).CoordID{1},sta_num)
                sat_sta_indx=ista;
            end
        end
        pb_match(day_count,isamp).sat_sta={all_sta(idt,sat_sta_indx).CoordID};
        pb_match(day_count,isamp).sat_lat={all_sta(idt,sat_sta_indx).Latitude};
        pb_match(day_count,isamp).sat_lon={all_sta(idt,sat_sta_indx).Longitude};
        pb_match(day_count,isamp).sat_date_time=all_sta(idt,sat_sta_indx).date_time;
        
        switch algthm
            case 'c2rcc'
                % C2RCC and C2X
                pb_match(day_count,isamp).sat_chl=all_sta(idt,sat_sta_indx).conc_chl;
                pb_match(day_count,isamp).sat_kd489=all_sta(idt,sat_sta_indx).kd489;
                pb_match(day_count,isamp).sat_kdz90=all_sta(idt,sat_sta_indx).kd_z90max;
                pb_match(day_count,isamp).sat_rhow_1=all_sta(idt,sat_sta_indx).rhow_1;
                pb_match(day_count,isamp).sat_rhow_2=all_sta(idt,sat_sta_indx).rhow_2;
                pb_match(day_count,isamp).sat_rhow_3=all_sta(idt,sat_sta_indx).rhow_3;
                pb_match(day_count,isamp).sat_rhow_4=all_sta(idt,sat_sta_indx).rhow_4;
                pb_match(day_count,isamp).sat_rhow_5=all_sta(idt,sat_sta_indx).rhow_5;
                pb_match(day_count,isamp).sat_tsm=all_sta(idt,sat_sta_indx).conc_tsm;
            case 'c2x'
                % C2RCC and C2X
                pb_match(day_count,isamp).sat_chl=all_sta(idt,sat_sta_indx).conc_chl;
                pb_match(day_count,isamp).sat_kd489=all_sta(idt,sat_sta_indx).kd489;
                pb_match(day_count,isamp).sat_kdz90=all_sta(idt,sat_sta_indx).kd_z90max;
                pb_match(day_count,isamp).sat_rhow_1=all_sta(idt,sat_sta_indx).rhow_1;
                pb_match(day_count,isamp).sat_rhow_2=all_sta(idt,sat_sta_indx).rhow_2;
                pb_match(day_count,isamp).sat_rhow_3=all_sta(idt,sat_sta_indx).rhow_3;
                pb_match(day_count,isamp).sat_rhow_4=all_sta(idt,sat_sta_indx).rhow_4;
                pb_match(day_count,isamp).sat_rhow_5=all_sta(idt,sat_sta_indx).rhow_5;
                pb_match(day_count,isamp).sat_tsm=all_sta(idt,sat_sta_indx).conc_tsm;
            case 'l2gen'
                % L2gen
                pb_match(day_count,isamp).sat_chl=all_sta(idt,sat_sta_indx).chlor_a;
                pb_match(day_count,isamp).sat_kd489=all_sta(idt,sat_sta_indx).Kd_490;
                pb_match(day_count,isamp).sat_rhow_1=all_sta(idt,sat_sta_indx).Rrs_443;
                pb_match(day_count,isamp).sat_rhow_2=all_sta(idt,sat_sta_indx).Rrs_482;
                pb_match(day_count,isamp).sat_rhow_3=all_sta(idt,sat_sta_indx).Rrs_561;
                pb_match(day_count,isamp).sat_rhow_4=all_sta(idt,sat_sta_indx).Rrs_655;
        end

        % In situ
        pb_match(day_count,isamp).insitu_date_time=pb_date_time(dt_indx(isamp));
        pb_match(day_count,isamp).insitu_totdepth=pb_totdep(dt_indx(isamp));
        pb_match(day_count,isamp).insitu_secchi=pb_secchi(dt_indx(isamp));
        pb_match(day_count,isamp).insitu_chl=pb_chl(dt_indx(isamp));
        pb_match(day_count,isamp).insitu_do=pb_do(dt_indx(isamp));
    end
end

save([outdir,'pb_insitu_sat_match_',num2str(hour_window),'h.mat'],"pb_match");

disp('Completed');