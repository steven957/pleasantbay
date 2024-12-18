% plot_annual_time_series_pb.m
% Syntax: plot_annual_time_series_pb
%
% Script loads satellite matchup data file and finds
% colocated data within a designated matchup window for
% time series over the study period (2013-2022)
%
% Inputs:
%    1) Directory location for satellite matchup file pb_all_sat.mat
%       file generated with read_matchup_pb.m
%
% Outputs:
%    1) Plots with annual time series for satellite observations
%   
% Other m-files required: None 
%
% MAT-files required: 
%    1) pb_all_sat.mat
%
% Author: Steven E. Lohrenz, Ph.D., biological oceanography
% School for Marine Science and Technology, University of Massachusetts Dartmouth
% email address: slohrenz@umassd.edu
% Website: http://www.umassd.edu/smast/
% Last revision: 20 Aug 2024

%% ------------- BEGIN CODE --------------%%Read matchup spreadsheet

clc
clearvars

%% Load data

% Load satellite matchup data
InDir1='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\ArcGIS\Projects\PleasantBay\Satellite_matchups\';
filename1='pb_all_sat.mat';
load([InDir1,filename1]);

% If using RF algorithm for Secchi, apply algorithm to rhow values
load([InDir1,'Matchup files\RF\','RF_chl_36h_Mdl.mat']);
load([InDir1,'Matchup files\RF\','RF_secchi_36h_Mdl.mat']);

year_range = 2013:2022;

season = 'summer';

% Set color order for plot
C = orderedcolors("gem12");
C2 = orderedcolors("glow12");
C = [C;C2];
colororder(C);

% Plot annual average time-series for each for each station for satellite and in situ data

year_chl_ts = nan(10,21);
year_kd_ts = nan(10,21);
year_RF_secchi = nan(10,21);
year_sta = cell(21,1);

% year_obs = nan(10,21);

for ista=1:21
    for iyr=1:length(year_range)
        switch season
            case 'summer'
                sta_indx = find(strcmp(pb_all_sat.sat_sta,num2str(ista)) & ...
                    isbetween(pb_all_sat.sat_date_time,datetime(year_range(iyr),6,1),...
                     datetime(year_range(iyr),8,eomday(year_range(iyr),8))));
            case 'spring'
                sta_indx = find(strcmp(pb_all_sat.sat_sta,num2str(ista)) & ...
                    isbetween(pb_all_sat.sat_date_time,datetime(year_range(iyr),3,1),...
                     datetime(year_range(iyr),9,eomday(year_range(iyr),5))));
            case 'fall'
                sta_indx = find(strcmp(pb_all_sat.sat_sta,num2str(ista)) & ...
                    isbetween(pb_all_sat.sat_date_time,datetime(year_range(iyr),9,1),...
                     datetime(year_range(iyr),11,eomday(year_range(iyr),11))));
            case 'winter'
                sta_indx = find(strcmp(pb_all_sat.sat_sta,num2str(ista)) & ...
                    isbetween(pb_all_sat.sat_date_time,datetime(year_range(iyr),12,1),...
                     datetime(year_range(iyr),2,eomday(year_range(iyr),2))));
        end                
        
        year_chl_ts(iyr,ista) = mean(pb_all_sat.sat_chl(sta_indx));
        year_kd_ts(iyr,ista) = mean(pb_all_sat.sat_kd489(sta_indx));
        
        % RF
        rhow_tab = table(pb_all_sat.sat_rhow_1(sta_indx),pb_all_sat.sat_rhow_2(sta_indx),...
            pb_all_sat.sat_rhow_3(sta_indx),pb_all_sat.sat_rhow_4(sta_indx),...
            pb_all_sat.sat_rhow_5(sta_indx),'VariableNames',["sat_rhow_440","sat_rhow_480",...
            "sat_rhow_560","sat_rhow_655","sat_rhow_865"]);

        RF_secchi = predict(Mdl,rhow_tab);
        RF_chl = predict(ChlMdl,rhow_tab);

        year_RF_secchi(iyr,ista) = mean(RF_secchi);
        year_RF_chl(iyr,ista) = mean(RF_chl);

        if strcmp(pb_all_sat.sat_sta(sta_indx),'2') | strcmp(pb_all_sat.sat_sta(sta_indx),'7') | ...
                strcmp(pb_all_sat.sat_sta(sta_indx),'18')
            year_chl_ts(iyr,ista)=nan;
            year_kd_ts(iyr,ista) = nan;
        end
    end
    year_sta(ista) = pb_all_sat.sat_sta(sta_indx(1));
end

%% Plot results for chl
% 
% Create a figure and adjust size to proportion of screen view
f1 = figure(1);
clf
scrsz = get(groot,'ScreenSize');
set(f1,'Position',[scrsz(4).*.3 scrsz(3).*.07 scrsz(3).*.45 scrsz(4).*.75])
set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

% For RF model, uncomment line below:
% year_chl_ts = year_RF_chl;

hp1=plot(year_range',year_chl_ts,'-','MarkerSize',8,'LineWidth',1.5);
box on
hold on
set(gca,'Fontsize',15,'Fontname','Arial','YLim',[0 18]);

for ista = 1:21
    htxt = text(year_range',year_chl_ts(:,ista),num2str(ista),'FontSize',12,'Color',C(ista,:),...
        'VerticalAlignment','bottom');
end

% yearlines = [2015,0;2015,10;nan,nan;2018,0;2018,10;nan,nan;2021,0;2021,10];
% hyearlines = plot(yearlines(:,1),yearlines(:,2),':k','LineWidth',1.5);

hyeararrows(1) = annotation('textarrow',[0.325 0.325], [0.97 0.87],'LineWidth',1.5,'Color','k');
% hyeararrows(2) = annotation('textarrow',[0.517 0.517], [0.91 0.81],'LineWidth',1.5,'Color','k');
hyeararrows(3) = annotation('textarrow',[0.711 0.711], [0.97 0.87],'LineWidth',1.5,'Color','k');

xlabel('Year','Fontsize',16,'Fontname','Arial','Fontweight','Bold')
ylabel('Landsat 8 C2RCC conc\_chl (mg m^{-3})', 'Fontsize',16,'Fontname','Arial','Fontweight','Bold');

% ylabel('Landsat 8 RF-Estimated Chl (mg m^{-3})', 'Fontsize',16,'Fontname','Arial','Fontweight','Bold');

% Prepare legend entries for each station
sta_char = num2str([1,3:6,8:17,19:21]');
hleg1 = legend(hp1([1,3:6,8:17,19:21]')',sta_char,'NumColumns',7,'Location','Southoutside');
% leg_pos = get(hleg1,'Position');
% leg_pos(1) = leg_pos(1) + 0.1;
% set(hleg1,"Position",leg_pos);

% print([InDir1,'c2rcc_chl_insitu_chl_annual_time_series_',season,'.tif'],'-r300','-dtiff');
% print([InDir1,'RF_chl_insitu_chl_annual_time_series_',season,'.tif'],'-r300','-dtiff');

% Evaluate correlation of conc_chl versus year
% First rearrange data 
chl_year_data = [repmat(year_range(2:end)',21,1),reshape(year_chl_ts(2:end,:),9.*21,1)];
chl_year_data = chl_year_data(~isnan(chl_year_data(:,2)),:);
[R_chl,P_chl] = corrcoef(chl_year_data);

%% Plot results for Secchi
% 
% Create a figure and adjust size to proportion of screen view
f3 = figure(2);
clf
scrsz = get(groot,'ScreenSize');
set(f3,'Position',[scrsz(4).*.3 scrsz(3).*.07 scrsz(3).*.45 scrsz(4).*.75])
set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

% Set color order for plot
C = orderedcolors("gem12");
C2 = orderedcolors("glow12");
C = [C;C2];
colororder(C);

% For RF model, uncomment line below:
%year_kd_ts = 1./year_RF_secchi;

% C2RCC Secchi
hp2=plot(year_range',1./year_kd_ts,'-','MarkerSize',8,'LineWidth',1.5); %,'Color',C);
hold on
set(gca,'Fontsize',15,'Fontname','Arial','YLim',[0.6 5.5]);

for ista = 1:21
    htxt = text(year_range',1./year_kd_ts(:,ista),num2str(ista),'FontSize',12,'Color',C(ista,:),...
        'VerticalAlignment','bottom');
end

hyeararrows(1) = annotation('textarrow',[0.325 0.325], [0.95 0.85],'LineWidth',1.5,'Color','k');
% hyeararrows(2) = annotation('textarrow',[0.517 0.517], [0.91 0.81],'LineWidth',1.5,'Color','k');
hyeararrows(3) = annotation('textarrow',[0.711 0.711], [0.95 0.85],'LineWidth',1.5,'Color','k');

xlabel('Year','Fontsize',16,'Fontname','Arial','Fontweight','Bold')

% RF Secchi
% ylabel('Landsat 8 RF-Estimated Secchi Depth (m)', 'Fontsize',16,'Fontname','Arial','Fontweight','Bold');

% C2RCC Secchi
ylabel('Landsat 8 C2RCC 1/kd489 (m)', 'Fontsize',16,'Fontname','Arial','Fontweight','Bold');

% Prepare legend entries for each station
sta_char = num2str([1,3:6,8:17,19:21]');
hleg1 = legend(hp2([1,3:6,8:17,19:21]')',sta_char,'NumColumns',7,'Location','Southoutside');
% leg_pos = get(hleg1,'Position');
% leg_pos(1) = leg_pos(1) + 0.1;
% set(hleg1,"Position",leg_pos);

% print([InDir1,'RF_secchi_insitu_secchi_annual_time_series_',season,'.tif'],'-r300','-dtiff');
% print([InDir1,'1overkd489_insitu_secchi_annual_time_series_',season,'.tif'],'-r300','-dtiff');

% Evaluate correlation of 1/kd489 versus year
% First rearrange data 
% inversekd_year_data = [repmat(year_range(2:end)',21,1),reshape(1./year_kd_ts(2:end,:),9.*21,1)];
% inversekd_year_data = inversekd_year_data(~isnan(inversekd_year_data(:,2)),:);
% [R_1overkd,P_1overkd] = corrcoef(inversekd_year_data);

% First rearrange data 
RF_secchi_year_data = [repmat(year_range(2:end)',21,1),reshape(year_RF_secchi(2:end,:),9.*21,1)];
RF_secchi_year_data = RF_secchi_year_data(~isnan(RF_secchi_year_data(:,2)),:);
[R_RF_secchi,P_RF_secchi] = corrcoef(RF_secchi_year_data);

%% Plot annual time series for precipitation
% 
% Load precipitation data from spreadsheet

precip_table = readtable([InDir1,'..\','plymouth_county_weather_data_2013-2024.xlsx']);

% Create a figure and adjust size to proportion of screen view
f5 = figure(3);
clf
scrsz = get(groot,'ScreenSize');
set(f5,'Position',[scrsz(4).*.3 scrsz(3).*.07 scrsz(3).*.45 scrsz(4).*.65])
set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

% Set color order for plot
C = orderedcolors("gem12");
C2 = orderedcolors("glow12");
C = [C;C2];

% Plot annual time series for each station for satellite and in situ data for each station and month
plt_dat = zeros(10,1);
for iyr = 1:length(year_range)
    switch season
        case 'spring'
            dat_indx = find((contains(precip_table.Month,'April') | contains(precip_table.Month,'March') | ...
                contains(precip_table.Month,'May')) & precip_table.Year==year_range(iyr));
        case 'summer'
            dat_indx = find((contains(precip_table.Month,'June') | contains(precip_table.Month,'July') | ...
                contains(precip_table.Month,'August') | contains(precip_table.Month,'September')) & ...
                precip_table.Year==year_range(iyr));
    end
    plt_dat(iyr) = mean(precip_table.Precipitation(dat_indx)).*2.54;  % Convert units from in to cm
end

hp5=bar(year_range,plt_dat,0.6);
% hp1=plot(year_range,plt_dat,'-o','MarkerSize',6,'LineWidth',1.5,'Color',C(1,:));
set(gca,'XLim',[2013.5 2022.5],'Fontsize',15,'Fontname','Arial');

xlabel('Year','Fontsize',16,'Fontname','Arial','Fontweight','Bold')
ylabel('Precipation, Plymouth County, MA (cm)', 'Fontsize',16,'Fontname','Arial','Fontweight','Bold');
% title('In Situ vs Satellite Chlorophyll-A', 'Fontsize',20,'Fontname','Arial','Fontweight','Bold');

hyeararrows(1) = annotation('textarrow',[0.3465 0.3465], [0.95 0.85],'LineWidth',1.5,'Color','k');
% hyeararrows(2) = annotation('textarrow',[0.517 0.517], [0.91 0.81],'LineWidth',1.5,'Color','k');
hyeararrows(3) = annotation('textarrow',[0.69 0.69], [0.95 0.85],'LineWidth',1.5,'Color','k');

print([InDir1,'plymouth_precipitation_annual_time_series.tif'],'-r300','-dtiff');

% Evaluate correlation of precipitation versus year
% First rearrange data 
precip_year_data = [year_range(2:end)',plt_dat(2:end)];
[R_precip,P_precip] = corrcoef(precip_year_data);


disp('Completed');
