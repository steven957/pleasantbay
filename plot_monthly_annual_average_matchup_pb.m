% plot_monthly_annual_average_matchup_pb.m
% Syntax: plot_monthly_annual_average_matchup_pb
%
% Script loads in_situ data and satellite matchup data files and finds
% colocated data within a designated matchup window for monthly
% climatologies over the study period (2013-2022)
%
% Inputs:
%    1) Directory location for satellite matchup all_station_match.mat file 
%       generated with matchup_PB_sat_insitu.mat and in situ matchup pb_match.mat file 
%       generated with matchup_PB_sat_insitu.mat
%
% Outputs:
%    1) Plots comparing monthly climatologies for satellite and in situ observations
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
% Last revision: 26 June 2024

%% ------------- BEGIN CODE --------------%%Read matchup spreadsheet
clc
clearvars


%% Load data

% Load satellite matchup data
InDir1='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\ArcGIS\Projects\PleasantBay\Satellite_matchups\';
filename1='all_station_matchup.mat';
load([InDir1,filename1]);

flds = fieldnames(all_sta);
[mimg,nsta] = size(all_sta);

% Load in situ data
InDir2='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\ArcGIS\Projects\PleasantBay\InSitu_Data\';
filename2='all_pb_dat_2010_2021.mat';
load([InDir2,filename2]);

% Create table for in situ data
all_pb_insitu = table(pb_date_time,pb_sta,pb_chl,pb_secchi,pb_do,pb_totdep);
all_pb_insitu.Properties.VariableNames = ["Date_Time","Station","Chl","Secchi","DO","Tot_Dep"];
[minsitu,~] = size(all_pb_insitu);

% Remove bad chlorophyll and secchi values
all_pb_insitu.Chl(all_pb_insitu.Chl==0.03,1)=nan;
all_pb_insitu.Secchi(all_pb_insitu.Secchi>=all_pb_insitu.Tot_Dep)=nan;

%% Compile satellite matchups for all dates for each station

for idy=1:mimg
    for ista=1:nsta
        if isempty(all_sta(idy,ista).CoordID)
            continue
        end

        % Find matching records for each day and station
        if ~exist('all_pb_sat','var')
            all_pb_sat=table(all_sta(idy,ista).date_time(1),all_sta(idy,ista).CoordID(1),mean(all_sta(idy,ista).conc_chl),...
            std(all_sta(idy,ista).conc_chl),mean(all_sta(idy,ista).kd489),...
            std(all_sta(idy,ista).kd489));
        else
            all_pb_sat=[all_pb_sat;{all_sta(idy,ista).date_time(1),all_sta(idy,ista).CoordID(1),...
            mean(all_sta(idy,ista).conc_chl),std(all_sta(idy,ista).conc_chl),mean(all_sta(idy,ista).kd489),...
            std(all_sta(idy,ista).kd489)}];
        end

        % Single pixel
        % all_chl=[all_chl;[pb_match(idy,ista).insitu_chl,pb_match(idy,ista).sat_chl(5)]];
        % all_kd489_kd489=[all_kd489_kd489;[pb_match(idy,ista).insitu_kd489,...
        %     1./pb_match(idy,ista).sat_kd489(5),pb_match(idy,ista).insitu_totdepth]];
        
    end
end

% Change variable names for table to human readable
all_pb_sat.Properties.VariableNames = ["Date_Time","Station","mean_conc_chl","std_conc_chl",...
    "mean_kd489","std_kd489"];

%% Find averages for for each station, month and year

year_range = 2013:2022;

sat_indx = cell(12,1);
obs_indx = cell(12,1);
for imnth=1:12
    for ista2=1:nsta
        % Find matching records for each month and station for a given
        %   year
        for iyr=1:length(year_range)
            sat_indx{imnth} = [sat_indx{imnth};find(isbetween(all_pb_sat.Date_Time,datetime(year_range(iyr),imnth,1),...
                datetime(year_range(iyr),imnth,eomday(year_range(iyr),imnth))) & ...
                strcmp(all_pb_sat.Station,num2str(ista2)))];
    
            obs_indx{imnth} = [obs_indx{imnth};find(isbetween(all_pb_insitu.Date_Time,datetime(year_range(iyr),imnth,1),...
                datetime(year_range(iyr),imnth,eomday(year_range(iyr),imnth))) & ...
                strcmp(all_pb_insitu.Station,['PBA-',num2str(ista2)]))];
                 % & ~strcmp(all_pb_insitu.Station,'PBA-1'))]; % & ~strcmp(all_pb_insitu.Station,'PBA-2'))];
        end

        % Create a table with monthly average data for each station and
        % assigning nan to variables without data for a given month
        if ~exist('all_pb_monthly_data','var') % Create table if does not already exist
            if isempty(obs_indx{imnth})  % Check for empty index array
                all_pb_monthly_data = table(imnth,ista2,mean(all_pb_sat.mean_conc_chl(sat_indx{imnth}),'omitmissing'),...
                    mean(all_pb_sat.mean_kd489(sat_indx{imnth}),'omitmissing'),...
                    std(all_pb_sat.mean_conc_chl(sat_indx{imnth}),'omitmissing'),...
                    std(all_pb_sat.mean_kd489(sat_indx{imnth}),'omitmissing'),nan,nan,nan,nan);
            else
                all_pb_monthly_data = table(imnth,ista2,mean(all_pb_sat.mean_conc_chl(sat_indx{imnth}),'omitmissing'),...
                    mean(all_pb_sat.mean_kd489(sat_indx{imnth}),'omitmissing'),...
                    std(all_pb_sat.mean_conc_chl(sat_indx{imnth}),'omitmissing'),...
                    std(all_pb_sat.mean_kd489(sat_indx{imnth}),'omitmissing'),...
                    mean(all_pb_insitu.Chl(obs_indx{imnth}),'omitmissing'),...
                    mean(all_pb_insitu.Secchi(obs_indx{imnth}),'omitmissing'),...
                    std(all_pb_insitu.Chl(obs_indx{imnth}),'omitmissing'),...
                    std(all_pb_insitu.Secchi(obs_indx{imnth}),'omitmissing'));
            end
        else
            if isempty(obs_indx{imnth}) % Check for empty index array
                all_pb_monthly_data = [all_pb_monthly_data;...
                    {imnth,ista2,mean(all_pb_sat.mean_conc_chl(sat_indx{imnth}),'omitmissing'),...
                    mean(all_pb_sat.mean_kd489(sat_indx{imnth}),'omitmissing'),...
                    std(all_pb_sat.mean_conc_chl(sat_indx{imnth}),'omitmissing'),...
                    std(all_pb_sat.mean_kd489(sat_indx{imnth}),'omitmissing'),nan,nan,nan,nan}];
            else
                all_pb_monthly_data = [all_pb_monthly_data;...
                    {imnth,ista2,mean(all_pb_sat.mean_conc_chl(sat_indx{imnth}),'omitmissing'),...
                    mean(all_pb_sat.mean_kd489(sat_indx{imnth}),'omitmissing'),...
                    std(all_pb_sat.mean_conc_chl(sat_indx{imnth}),'omitmissing'),...
                    std(all_pb_sat.mean_kd489(sat_indx{imnth}),'omitmissing'),...
                    mean(all_pb_insitu.Chl(obs_indx{imnth}),'omitmissing'),...
                    mean(all_pb_insitu.Secchi(obs_indx{imnth}),'omitmissing'),...
                    std(all_pb_insitu.Chl(obs_indx{imnth}),'omitmissing'),...
                    std(all_pb_insitu.Secchi(obs_indx{imnth}),'omitmissing')}];
            end
        end
    end
end

all_pb_monthly_data.Properties.VariableNames = ["Month","Station","conc_chl_mean","kd489_mean",...
    "conc_chl_std","kd489_std","insitu_chl_mean","insitu_secchi_mean","insitu_chl_std","insitu_secchi_std"];

save([InDir1,'all_pb_monthly_data_C2RCC.mat'],'all_pb_monthly_data');

%% Plot results for chl
% 
% Create a figure and adjust size to proportion of screen view
f1 = figure(1);
clf
scrsz = get(groot,'ScreenSize');
set(f1,'Position',[scrsz(4).*.3 scrsz(3).*.07 scrsz(3).*.45 scrsz(4).*.75])
set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

% Set color order for plot
C = orderedcolors("gem12");
C2 = orderedcolors("glow12");
C = [C;C2];

% Plot each station for satellite and in situ data for each month
for ista=[1,3:6,8:17,19:21]
    mnth_indx = find(all_pb_monthly_data.Station==ista);
    % hp1(ista)=errorbar(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.conc_chl_mean(mnth_indx),...
    %     all_pb_monthly_data.conc_chl_std(mnth_indx),'-x','MarkerSize',5,'LineWidth',1.5);
    hp1(ista)=plot(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.conc_chl_mean(mnth_indx),...
        '-','MarkerSize',8,'LineWidth',1.25,'Color',C(ista,:));
    if ista==1
        hold on
    end

    % hp1=errorbar(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.insitu_chl_mean(mnth_indx),...
    %     all_pb_monthly_data.insitu_chl_std(mnth_indx),':o','MarkerSize',5,'LineWidth',1.5);
    hp2(ista)=plot(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.insitu_chl_mean(mnth_indx),...
        ':o','MarkerSize',8,'LineWidth',1.5,'Color',C(ista,:));

    text(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.conc_chl_mean(mnth_indx),num2str(ista),'FontSize',12,'Color',C(ista,:),...
        'VerticalAlignment','bottom');
end
set(gca,'YLim',[0,8],'XLim',[0.5,12.5],'Fontsize',15,'Fontname','Arial');

% str1 = {['N=',num2str(chl_num_obs), '; r^{2}=',num2str(mean(chl_rsquare),3),'; {\itp} = ',num2str(mean(chl_prob),3)]};
% str2 = {['RMSE = ',num2str(sqrt(mean(chl_root_mean_sq_error)),3)]};
% text(1,29,str1);
% text(1,27,str2);

xlabel('Month','Fontsize',16,'Fontname','Arial','Fontweight','Bold')
ylabel('Landsat 8 C2RCC conc\_chl or In Situ Chl (mg m^{-3})', 'Fontsize',16,'Fontname','Arial','Fontweight','Bold');
% title('In Situ vs Satellite Chlorophyll-A', 'Fontsize',20,'Fontname','Arial','Fontweight','Bold');

% Prepare legend entries for each station
sta_char = num2str([1,3:6,8:17,19:21]');
hleg1 = legend(hp1([1,3:6,8:17,19:21])',sta_char,'NumColumns',7,'Location','Southoutside');
% leg_pos = get(hleg1,'Position');
% leg_pos(1) = leg_pos(1) + 0.1;
% set(hleg1,"Position",leg_pos);

print([InDir1,'c2rcc_chl_insitu_chl_monthly_time_series.tif'],'-r300','-dtiff');

% Create plot of satellite versus in situ monthly averages
figure(2);
clf
% hp3=errorbar(all_pb_monthly_data.insitu_chl_mean,all_pb_monthly_data.conc_chl_mean,...
%     all_pb_monthly_data.conc_chl_std,'ko', 'MarkerSize', 5);
plt_indx = find(all_pb_monthly_data.Station~=2 | all_pb_monthly_data.Station~=7 | all_pb_monthly_data.Station~=18);
hp3=plot(all_pb_monthly_data.insitu_chl_mean,all_pb_monthly_data.conc_chl_mean,'ko', 'MarkerSize', 5);
hold on 
plt_indx2 = find(all_pb_monthly_data.Station==1);
hp3=plot(all_pb_monthly_data.insitu_chl_mean(plt_indx2),all_pb_monthly_data.conc_chl_mean(plt_indx2),...
    'bo', 'MarkerSize',10,'LineWidth',1.5);
axpos = get(gca,'Position');
axpos(2)=axpos(2)+0.04;
set(gca,'Position',axpos,'XLim',[0 6],'YLim',[0 6],'Fontsize',14,'Fontname','Arial');  %,'XScale','log','YScale','log');
plot([0,6],[0,6],'k:');
xlabel('In Situ Chl (mg m^{-3})','Fontsize',15,'Fontname','Arial','Fontweight','Bold')
ylabel('Landsat 8 C2RCC conc\_chl (mg m^{-3})', 'Fontsize',15,'Fontname','Arial','Fontweight','Bold');

% Exclude Station 1 from fit
fit_indx = setdiff(plt_indx,plt_indx2); % plt_indx;  
lm2 = fitlm(all_pb_monthly_data.insitu_chl_mean(fit_indx),all_pb_monthly_data.conc_chl_mean(fit_indx));
chl_rsquare=lm2.Rsquared.Adjusted;
chl_root_mean_sq_error=lm2.RMSE;
chl_num_obs=lm2.NumObservations;
chl_prob=lm2.Coefficients.pValue(2);
sat_chl_est = predict(lm2,(0:20)');
% plot(0:20,sat_chl_est,'k-','Linewidth',2);

str1 = {['N=',num2str(chl_num_obs), '; r^{2}=',num2str(chl_rsquare,3),'; {\itp} = ',num2str(chl_prob,3)]};
str2 = {['RMSE = ',num2str(chl_root_mean_sq_error,3)]};
text(0.225,5.80,str1);
text(0.225,5.520,str2);

print([InDir1,'insitu_chl_vs_sat_conc_chl_monthly.tif'],'-r300','-dtiff');

%% Plot results for Secchi
% 
% Create a figure and adjust size to proportion of screen view
f3 = figure(3);
clf
scrsz = get(groot,'ScreenSize');
set(f3,'Position',[scrsz(4).*.3 scrsz(3).*.07 scrsz(3).*.45 scrsz(4).*.75])
set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

% Set color order for plot
C = orderedcolors("gem12");
C2 = orderedcolors("glow12");
C = [C;C2];

% Plot each station for satellite and in situ data for each month
for ista=[1,3:6,8:17,19:21]
    mnth_indx = find(all_pb_monthly_data.Station==ista);
    % hp1(ista)=errorbar(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.conc_chl_mean(mnth_indx),...
    %     all_pb_monthly_data.conc_chl_std(mnth_indx),'-x','MarkerSize',5,'LineWidth',1.5);
    hp4(ista)=plot(all_pb_monthly_data.Month(mnth_indx),1./all_pb_monthly_data.kd489_mean(mnth_indx),...
        '-','MarkerSize',8,'LineWidth',1.25,'Color',C(ista,:));
    if ista==1
        hold on
    end
    % hp1=errorbar(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.insitu_chl_mean(mnth_indx),...
    %     all_pb_monthly_data.insitu_chl_std(mnth_indx),':o','MarkerSize',5,'LineWidth',1.5);
    hp5(ista)=plot(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.insitu_secchi_mean(mnth_indx),...
        ':o','MarkerSize',8,'LineWidth',1.5,'Color',C(ista,:));
    
    text(all_pb_monthly_data.Month(mnth_indx),1./all_pb_monthly_data.kd489_mean(mnth_indx),num2str(ista),'FontSize',12,'Color',C(ista,:),...
        'VerticalAlignment','bottom');
end
set(gca,'YLim',[0.5,3],'XLim',[0.5,12.5],'Fontsize',15,'Fontname','Arial');

% str1 = {['N=',num2str(chl_num_obs), '; r^{2}=',num2str(mean(chl_rsquare),3),'; {\itp} = ',num2str(mean(chl_prob),3)]};
% str2 = {['RMSE = ',num2str(sqrt(mean(chl_root_mean_sq_error)),3)]};
% text(1,29,str1);
% text(1,27,str2);

xlabel('Month','Fontsize',16,'Fontname','Arial','Fontweight','Bold')
ylabel('Landsat 8 C2RCC 1/kd489 or Secchi Depth (m))', 'Fontsize',16,'Fontname','Arial','Fontweight','Bold');
% title('In Situ vs Satellite Chlorophyll-A', 'Fontsize',20,'Fontname','Arial','Fontweight','Bold');

% Prepare legend entries for each station
sta_char = num2str([1,3:6,8:17,19:21]');
hleg1 = legend(hp4([1,3:6,8:17,19:21])',sta_char,'NumColumns',7,'Location','Southoutside');
% leg_pos = get(hleg1,'Position');
% leg_pos(1) = leg_pos(1) + 0.1;
% set(hleg1,"Position",leg_pos);

print([InDir1,'1overkd489_insitu_secchi_monthly_time_series.tif'],'-r300','-dtiff');

% Create plot of satellite versus in situ monthly averages
figure(4);
clf
% hp3=errorbar(all_pb_monthly_data.insitu_secchi_mean,all_pb_monthly_data.kd489_mean,...
%     all_pb_monthly_data.kd489_std,'ko', 'MarkerSize', 5);
plt_indx = find(all_pb_monthly_data.Station~=2 | all_pb_monthly_data.Station~=7 | all_pb_monthly_data.Station~=18);
hp3=plot(all_pb_monthly_data.insitu_secchi_mean,1./all_pb_monthly_data.kd489_mean,'ko', 'MarkerSize', 5);
hold on 
plt_indx2 = find(all_pb_monthly_data.Station==1);
hp6=plot(all_pb_monthly_data.insitu_secchi_mean(plt_indx2),1./all_pb_monthly_data.kd489_mean(plt_indx2),...
    'bo', 'MarkerSize',10,'LineWidth',1.5);
axpos = get(gca,'Position');
axpos(2)=axpos(2)+0.04;
set(gca,'Position',axpos,'XLim',[1 3.5],'YLim',[1 3.5],'Fontsize',14,'Fontname','Arial');  %,'XScale','log','YScale','log');
plot([1 3.5],[1 3.5],'k:');
xlabel('Secchi Depth (m))','Fontsize',15,'Fontname','Arial','Fontweight','Bold')
ylabel('Landsat 8 C2RCC 1/kd489 (m)', 'Fontsize',15,'Fontname','Arial','Fontweight','Bold');
% title('Secchi vs Kd489', 'Fontsize',20,'Fontname','Arial','Fontweight','Bold');
lm2 = fitlm(all_pb_monthly_data.insitu_secchi_mean(plt_indx),1./all_pb_monthly_data.kd489_mean(plt_indx));
secchi_rsquare=lm2.Rsquared.Adjusted;
secchi_root_mean_sq_error=lm2.RMSE;
secchi_num_obs=lm2.NumObservations;
secchi_prob=lm2.Coefficients.pValue(2);
sat_secchi_est = predict(lm2,(0:20)');
% plot(0:20,sat_secchi_est,'k-','Linewidth',2);

str3 = {['N=',num2str(secchi_num_obs), '; r^{2}=',num2str(secchi_rsquare,3),'; {\itp} = ',num2str(secchi_prob,3)]};
str4 = {['RMSE = ',num2str(secchi_root_mean_sq_error,3)]};
text(1.05,3.425,str3);
text(1.05,3.3,str4);

print([InDir1,'insitu_secchi_vs_sat_kd489_monthly.tif'],'-r300','-dtiff');

%% Plot annual time series for each month
% 
% Create a figure and adjust size to proportion of screen view
%{
f5 = figure(5);
clf
scrsz = get(groot,'ScreenSize');
set(f5,'Position',[scrsz(4).*.3 scrsz(3).*.07 scrsz(3).*.45 scrsz(4).*.75])
set(0,'DefaultFigureVisible','on');  %Suppresses figure visibility during processing - set to on if desired

% Set color order for plot
C = orderedcolors("gem12");
C2 = orderedcolors("glow12");
C = [C;C2];

% Plot each station for satellite and in situ data for each month
for ista=1:nsta
    for iyr=1:length(year_range)
        yr_indx = find(all_pb_monthly_data.Station==ista);
        % hp1(ista)=errorbar(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.conc_chl_mean(mnth_indx),...
        %     all_pb_monthly_data.conc_chl_std(mnth_indx),'-x','MarkerSize',5,'LineWidth',1.5);
        hp1(ista)=plot(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.conc_chl_mean(mnth_indx),...
            '-x','MarkerSize',8,'LineWidth',1.5,'Color',C(ista,:));
        if ista==1
            hold on
        end
        % hp1=errorbar(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.insitu_chl_mean(mnth_indx),...
        %     all_pb_monthly_data.insitu_chl_std(mnth_indx),':o','MarkerSize',5,'LineWidth',1.5);
        hp2(ista)=plot(all_pb_monthly_data.Month(mnth_indx),all_pb_monthly_data.insitu_chl_mean(mnth_indx),...
            ':o','MarkerSize',8,'LineWidth',1.5,'Color',C(ista,:));
    end
end
set(gca,'YLim',[0,8],'XLim',[0.5,12.5],'Fontsize',15,'Fontname','Arial');

% str1 = {['N=',num2str(chl_num_obs), '; r^{2}=',num2str(mean(chl_rsquare),3),'; {\itp} = ',num2str(mean(chl_prob),3)]};
% str2 = {['RMSE = ',num2str(sqrt(mean(chl_root_mean_sq_error)),3)]};
% text(1,29,str1);
% text(1,27,str2);

xlabel('Month','Fontsize',16,'Fontname','Arial','Fontweight','Bold')
ylabel('Landsat 8 C2RCC conc\_chl or In Situ Chl (mg m^{-3})', 'Fontsize',16,'Fontname','Arial','Fontweight','Bold');
% title('In Situ vs Satellite Chlorophyll-A', 'Fontsize',20,'Fontname','Arial','Fontweight','Bold');

% Prepare legend entries for each station
sta_char = num2str((1:21)');
hleg1 = legend(hp1',sta_char,'NumColumns',7,'Location','Southoutside');
% leg_pos = get(hleg1,'Position');
% leg_pos(1) = leg_pos(1) + 0.1;
% set(hleg1,"Position",leg_pos);

print([InDir1,'..\','c2rcc_chl_insitu_chl_monthly_time_series.tif'],'-r300','-dtiff');
%}

disp('Completed');
