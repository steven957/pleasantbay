% plot_matchup_pb.m
% Syntax: plot_matchup_pb
%
% Script loads in_situ data and satellite matchup data files and finds
% colocated data within a designated matchup window
%
% Inputs:
%    1) Directory location for satellite and in situ matchup .mat file 
%       generated with matchup_PB_sat_insitu.mat
%
% Outputs:
%    1) Plots comparing satellite and in situ observations with options for
%    different algorithms
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

%% ------------- BEGIN CODE --------------

% Read matchup spreadsheet
clc
clearvars

algthm = 'c2rcc';

switch algthm
    case 'c2rcc'
        %C2RCC
        InDir='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\ArcGIS\Projects\PleasantBay\Satellite_matchups\';
    case 'c2x'
        % C2X
        InDir='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\ArcGIS\Projects\PleasantBay\Satellite_matchups\Matchup files\C2X\';
    case 'l2gen'
        % L2gen
        InDir='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\ArcGIS\Projects\PleasantBay\Satellite_matchups\Matchup files\l2gen\';
end

%Load data

matchup_window = '6h';  % '6h';

filename=['pb_insitu_sat_match_',matchup_window,'.mat'];

load([InDir,filename]);

flds = fieldnames(pb_match);
[mimg,nsta] = size(pb_match);

all_chl = [];
all_secchi_kd489 = [];
for idy=1:mimg
    for ista=1:nsta
        if isempty(pb_match(idy,ista).insitu_sta) || strcmp(pb_match(idy,ista).insitu_sta,'2')
            continue
        end
        all_chl=[all_chl;[pb_match(idy,ista).insitu_chl,mean(pb_match(idy,ista).sat_chl),...
            std(pb_match(idy,ista).sat_chl)]];
        switch algthm
            case 'c2rcc'
                all_secchi_kd489=[all_secchi_kd489;[pb_match(idy,ista).insitu_secchi,...
                mean(1./pb_match(idy,ista).sat_kd489),std(1./pb_match(idy,ista).sat_kd489),...
                pb_match(idy,ista).insitu_totdepth,mean(pb_match(idy,ista).sat_kdz90),...
                std(pb_match(idy,ista).sat_kdz90)]];
            case 'c2x'
                all_secchi_kd489=[all_secchi_kd489;[pb_match(idy,ista).insitu_secchi,...
                mean(1./pb_match(idy,ista).sat_kd489),std(1./pb_match(idy,ista).sat_kd489),...
                pb_match(idy,ista).insitu_totdepth,mean(pb_match(idy,ista).sat_kdz90),...
                std(pb_match(idy,ista).sat_kdz90)]];
            case 'l2gen'
                all_secchi_kd489=[all_secchi_kd489;[pb_match(idy,ista).insitu_secchi,...
                mean(1./pb_match(idy,ista).sat_kd489),std(1./pb_match(idy,ista).sat_kd489),...
                pb_match(idy,ista).insitu_totdepth]];
        end
                      
        % Single pixel
        % all_chl=[all_chl;[pb_match(idy,ista).insitu_chl,pb_match(idy,ista).sat_chl(5)]];
        % all_secchi_kd489=[all_secchi_kd489;[pb_match(idy,ista).insitu_secchi,...
        %     1./pb_match(idy,ista).sat_kd489(5),pb_match(idy,ista).insitu_totdepth]];
        
        % Remove bad chlorophyll and secchi values (bad chlorophylls were a
        %    batch of bad analyses with anomalously low values; bad Secchis
        %    were stations where inverse kd489 depth was greater than or equal to total
        %    depth)
        all_chl(all_chl(:,1)<=0.03 | all_chl(:,2)>100 | all_chl(:,2)<=0,:)=nan;
        all_secchi_kd489(all_secchi_kd489(:,1)>=all_secchi_kd489(:,4),1)=nan;
        % all_secchi_kd489(all_secchi_kd489(:,1)>=all_secchi_kd489(:,3),1)=nan;  % Single pixel
    end
end

figure(1);
clf
hp1=errorbar(all_chl(:,1),all_chl(:,2),all_chl(:,3),'ko','MarkerSize',5);
% hp1=plot(all_chl(:,1),all_chl(:,2),'ko', 'MarkerSize', 5);  % Single pixel

hold on
plot([.1,30],[.1,30],'k:');
lm1 = fitlm(all_chl(:,1),all_chl(:,2));
chl_rsquare=lm1.Rsquared.Adjusted;
chl_root_mean_sq_error=lm1.RMSE;
chl_num_obs=lm1.NumObservations;
chl_prob=lm1.Coefficients.pValue(2);
sat_chl_est = predict(lm1,(0:20)');

str1 = {['N=',num2str(chl_num_obs), '; r^{2}=',num2str(mean(chl_rsquare),3),'; {\itp} = ',num2str(mean(chl_prob),3)]};
str2 = {['RMSE = ',num2str(chl_root_mean_sq_error,3)]};

set(gca,'FontSize',14);

switch algthm
    case 'c2rcc'
        % C2RCC
        text(1,29,str1,'FontSize',11);
        text(1,27,str2,'FontSize',11);
        xlabel('In Situ Chla (mg m^{-3})','Fontsize',14,'Fontname','Arial','Fontweight','Bold')
        ylabel('Landsat 8 C2RCC conc\_chl (mg m^{-3})', 'Fontsize',14,'Fontname','Arial','Fontweight','Bold');
    case 'c2x'
        text(1,29,str1,'FontSize',11);
        text(1,27,str2,'FontSize',11);
        xlabel('In Situ Chla (mg m^{-3})','Fontsize',14,'Fontname','Arial','Fontweight','Bold')
        ylabel('Landsat 8 C2X conc\_chl (mg m^{-3})', 'Fontsize',14,'Fontname','Arial','Fontweight','Bold');
    case 'l2gen'        
        % L2gen
        text(1,80,str1,'FontSize',11);
        text(1,62,str2,'FontSize',11);
        set(gca,'YScale','log','XScale','log','YLim',[1,110]);
        xlabel('In Situ Chla (mg m^{-3})','Fontsize',14,'Fontname','Arial','Fontweight','Bold')
        ylabel('Landsat 8 L2gen chlor\_a (mg m^{-3})', 'Fontsize',14,'Fontname','Arial','Fontweight','Bold');
end

print([InDir,algthm,'_chl_vs_insitu_',matchup_window,'.tif'],'-r300','-dtiff');

figure(2);
clf
hp2=errorbar(all_secchi_kd489(:,1),all_secchi_kd489(:,2),all_secchi_kd489(:,3),'ko', 'MarkerSize', 5);
% hp2=plot(all_secchi_kd489(:,1),all_secchi_kd489(:,2),'ko', 'MarkerSize', 5); % Single pixel
hold on 
axpos = get(gca,'Position');
axpos(2)=axpos(2)+0.04;
set(gca,'Position',axpos); %,'XScale','log','YScale','log');
set(gca,'FontSize',14,'YLim',[0 6]); %,'XLim',[0 4.5]);
plot([0,5],[0,5],'k:');
xlabel('Secchi depth (m^{-1})','Fontsize',14,'Fontname','Arial','Fontweight','Bold')

switch algthm
    case 'c2rcc'
        ylabel('Landsat 8 C2RCC 1/kd489 (m)', 'Fontsize',14,'Fontname','Arial','Fontweight','Bold');
    case 'c2x'
        ylabel('Landsat 8 C2X 1/kd489 (m)', 'Fontsize',14,'Fontname','Arial','Fontweight','Bold');
    case 'l2gen'
        ylabel('Landsat 8 L2gen 1/Kd_490 (m)', 'Fontsize',14,'Fontname','Arial','Fontweight','Bold');
end
lm2 = fitlm(all_secchi_kd489(~isinf(all_secchi_kd489(:,2)),1),all_secchi_kd489(~isinf(all_secchi_kd489(:,2)),2));
secchi_rsquare=lm2.Rsquared.Adjusted;
secchi_root_mean_sq_error=lm2.RMSE;
secchi_num_obs=lm2.NumObservations;
secchi_prob=lm2.Coefficients.pValue(2);
sat_secchi_est = predict(lm2,(0:20)');
% plot(0:20,sat_secchi_est,'k-','Linewidth',2);

str1 = {['N=',num2str(secchi_num_obs), '; r^{2}=',num2str(mean(secchi_rsquare),3),'; {\itp} = ',num2str(mean(secchi_prob),3)]};
str2 = {['RMSE = ',num2str(secchi_root_mean_sq_error,3)]};
text(0.225,5.80,str1,'FontSize',11);
text(0.225,5.49,str2,'FontSize',11);

print([InDir,'secchi_vs_1over',algthm,'kd489_',matchup_window,'.tif'],'-r300','-dtiff');

switch algthm
    case 'l2gen'
        return  % C2RCC and C2X only
end

figure(3);
clf
hp2=errorbar(all_secchi_kd489(:,1),all_secchi_kd489(:,5),all_secchi_kd489(:,6),'ko', 'MarkerSize', 5);
% hp2=plot(all_secchi_kd489(:,1),all_secchi_kd489(:,2),'ko', 'MarkerSize', 5); % Single pixel
hold on 
axpos = get(gca,'Position');
axpos(2)=axpos(2)+0.04;
set(gca,'Position',axpos); %,'XScale','log','YScale','log');
set(gca,'FontSize',14,'YLim',[0 7]); %,'XLim',[0 4.5]);
plot([0,5],[0,5],'k:');
xlabel('Secchi depth (m^{-1})','Fontsize',14,'Fontname','Arial','Fontweight','Bold')

lm2 = fitlm(all_secchi_kd489(:,1),all_secchi_kd489(:,5));
secchi_kdz_rsquare=lm2.Rsquared.Adjusted;
secchi_kdz_root_mean_sq_error=lm2.RMSE;
secchi_kdz_num_obs=lm2.NumObservations;
secchi_kdz_prob=lm2.Coefficients.pValue(2);
sat_secchi_est = predict(lm2,(0:20)');
% plot(0:20,sat_secchi_est,'k-','Linewidth',2);

str1 = {['N=',num2str(secchi_num_obs), '; r^{2}=',num2str(mean(secchi_rsquare),3),'; {\itp} = ',num2str(mean(secchi_prob),3)]};
str2 = {['RMSE = ',num2str(secchi_root_mean_sq_error,3)]};
text(2.5,0.9,str1,'FontSize',11);
text(2.5,0.5,str2,'FontSize',11);

print([InDir,'secchi_vs_kd_z90max_',matchup_window,'.tif'],'-r300','-dtiff');

disp('Completed');
