% Pleasant_Bay_RndmFrst_matchup.m
% Syntax:  Pleasant_Bay_RndmFrst_matchup
%
% Script reads prediction data from spreadsheet and applies random forest
% regression to estimate water quality variables
%
% Inputs:
%    1) Directory locations for matchup files
%
% Outputs:
%    1) Plot of model regression fit
%    2) Plot of predictor importance
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
% Last revision: 23 Dec 2024

%% ------------- BEGIN CODE --------------%

clearvars
close all
clc

InDir='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\ArcGIS\Projects\PleasantBay\Satellite_matchups\';
OutDir='C:\Users\slohrenz\OneDrive - UMASS Dartmouth\Documents\ArcGIS\Projects\PleasantBay\';

%Load data
matchup_window = '36h';  % '36h';
filename=['pb_insitu_sat_match_',matchup_window,'.mat'];

rand('seed',0); % set random seed so results are reproducible

% Variable flag
plt_var = 'chl';  %'secchi'; %'chl';  %  'DO';

% Save plot flag
plt_save = 0;

% Save model flag
mdl_flag = 1;

load([InDir,filename]);

flds = fieldnames(pb_match);
[mimg,nsta] = size(pb_match);

count = 0;
for idy=1:mimg
    for ista=1:nsta
        if isempty(pb_match(idy,ista).insitu_sta) || strcmp(pb_match(idy,ista).insitu_sta,'2')
            continue
        end
        count = count + 1;
        insitu_sta(count,1) = pb_match(idy,ista).insitu_sta;
        sat_date_time(count,1) = pb_match(idy,ista).sat_date_time(1);

        insitu_totdepth(count,1) = pb_match(idy,ista).insitu_totdepth;
        insitu_chl(count,1) = pb_match(idy,ista).insitu_chl;
        insitu_secchi(count,1) = pb_match(idy,ista).insitu_secchi;
        insitu_do(count,1) =  pb_match(idy,ista).insitu_do;

        sat_chl(count,1) = mean(pb_match(idy,ista).sat_chl);
        inverse_kd489(count,1) = mean(1./pb_match(idy,ista).sat_kd489);
        sat_rhow_440(count,1) = mean(pb_match(idy,ista).sat_rhow_1);
        sat_rhow_480(count,1) = mean(pb_match(idy,ista).sat_rhow_2);
        sat_rhow_560(count,1) = mean(pb_match(idy,ista).sat_rhow_3);
        sat_rhow_655(count,1) = mean(pb_match(idy,ista).sat_rhow_4);
        sat_rhow_865(count,1) = mean(pb_match(idy,ista).sat_rhow_5);
        sat_tsm(count,1) = mean(pb_match(idy,ista).sat_tsm);
    end

        % Single pixel 440	480	560	655	865

        % all_chl=[all_chl;[pb_match(idy,ista).insitu_chl,pb_match(idy,ista).sat_chl(5)]];
        % all_secchi_kd489=[all_secchi_kd489;[pb_match(idy,ista).insitu_secchi,...
        %     1./pb_match(idy,ista).sat_kd489(5),pb_match(idy,ista).insitu_totdepth]];
end

all_data=table(insitu_sta,sat_date_time,insitu_totdepth,insitu_chl,...
    insitu_secchi,insitu_do,sat_chl,sat_rhow_440,sat_rhow_480,...
    sat_rhow_560,sat_rhow_655,sat_rhow_865,sat_tsm);

% Remove bad chlorophyll and secchi values
% all_data(all_data(:,1)==0.03,1)=nan;
% all_secchi_kd489(all_secchi_kd489(:,1)>=all_secchi_kd489(:,4),1)=nan;
% all_secchi_kd489(all_secchi_kd489(:,1)>=all_secchi_kd489(:,3),1)=nan;  % Single pixel

Station=all_data.insitu_sta;
R440=all_data.sat_rhow_440;
R480=all_data.sat_rhow_480;
R560=all_data.sat_rhow_560;
R655=all_data.sat_rhow_655;
R865=all_data.sat_rhow_865;


switch plt_var
    case 'chl'
        % For chlorophyll:
        IS=all_data.insitu_chl;
        good_indx = find(IS~=0.03);
    case 'secchi'
        % For secchi:
        IS=all_data.insitu_secchi;
        good_indx = find(IS<=all_data.insitu_totdepth);
    case 'DO'
        % For dissolved oxygen:
        IS=all_data.insitu_do;
        good_indx = find(IS>0);
end

tbl=all_data(good_indx,8:12);

Y=IS(good_indx); % Response variable dataset
tbl.Y=Y; % Add to evaluations dataset
n=length(Y);

nfolds=5;
hpartition=cvpartition(n,'Kfold',nfolds);  %Non-stratified partition (appropriate for regression problem)

% Develop five sets of training and test data
preY=[];
insituY=[];
trainError=zeros(5,1);
cvtrainError=zeros(5,1);
newError=zeros(5,1);
rsquare=zeros(5,1);
prob=zeros(5,1);
pred_imp = zeros(5,5);
for irep=1:nfolds
    idxTrain = training(hpartition,irep);
    tblTrain = tbl(idxTrain,:);
    idxNew = test(hpartition,irep);
    tblNew = tbl(idxNew,:);
    
    Mdl = fitrtree(tblTrain,'Y');  % Train random forest model with training data
    trainError(irep) = resubLoss(Mdl);  % Mean squared error for resubsitution
    
    cvMdl = crossval(Mdl); % Performs stratified 10-fold cross-validation
    cvtrainError(irep) = kfoldLoss(cvMdl); % Mean squared error for misclassification
    
    newError(irep) = loss(Mdl,tblNew,'Y');  % Not sure how this is used

    preY=[preY;predict(Mdl,tblNew(:,1:5))];
    insituY=[insituY;tblNew.Y];

    [corr,pro]=corrcoef(insituY,preY,'Rows','complete');
    rsquare(irep)=corr(2,1).^2;
    prob(irep)=pro(2,1);

    pred_imp(irep,:) = predictorImportance(Mdl);
end

switch plt_var
    case 'chl'
        if mdl_flag
            % Save RF model output
            ChlMdl = Mdl;
            save([InDir,'Matchup files\RF\','RF_chl_36h_Mdl.mat'],"ChlMdl");
        end

        %% Plot chlorophyll results
        
        figure(1);
        clf
        plot(insituY,preY,'ko');
        hold on
        xlabel('RF Predicted Chlorophyll (mg/m^3)','Fontsize',10,'Fontname','Arial','Fontweight','Bold')
        ylabel('In-Situ Chlorophyll (mg/m^3)', 'Fontsize',10,'Fontname','Arial','Fontweight','Bold');
        set(gca,'Fontsize',14,'Xlim',[0.4,60],'Ylim',[0.4,60],'XScale','log','YScale','log',...
            'XTick',[0.5,1.0,2.0,5.0,10,20,50],'YTick',[0.5,1.0,2.0,5.0,10,20,50]);
        box on
        plot(gca,[0.6,30],[0.6,30],':k','Linewidth',1.5);
        str1 = {['N=',num2str(n), '; r^{2}=',num2str(mean(rsquare),3),'; {\itp} = ',num2str(mean(prob),3)]};
        str2 = {['Substitution Error = ',num2str(sqrt(mean(trainError)),3), '; Classification Error =', num2str(sqrt(mean(cvtrainError)),3)]};
        text(0.5,50,str1);
        text(0.5,38,str2);
        
        if plt_save
            print(gcf,[OutDir,'RF_regression_non-stratified_partition_Chl_',matchup_window,'.tif'],'-dtiff','-r600');
        end
        
        figure(2);
        clf
        hp2 = bar({'rhow\_440','rhow\_480','rhow\_560','rhow\_655','rhow\_865'},mean(pred_imp,2),0.6);
        hold on
        xlabel('Predictors','Fontsize',10,'Fontname','Arial','Fontweight','Bold')
        ylabel('Predictor Importance', 'Fontsize',10,'Fontname','Arial','Fontweight','Bold');
        set(gca,'Fontsize',14);
        box on
        
        if plt_save
            print(gcf,[OutDir,'RF_regression_predictor_importance_Chl_',matchup_window,'.tif'],'-dtiff','-r600');
        end

    case 'secchi'
        if mdl_flag
            % Save RF model output
            secchiMdl = Mdl;
            save([InDir,'Matchup files\RF\','RF_secchi_36h_Mdl.mat'],"secchiMdl");
        end
        %% Plot Secchi results

        figure(1);
        clf
        hp1 = plot(insituY,preY,'ko');
        hold on
        xlabel('RF Predicted Secchi Depth (m)','Fontsize',10,'Fontname','Arial','Fontweight','Bold')
        ylabel('In-Situ Secchi Depth (m)', 'Fontsize',10,'Fontname','Arial','Fontweight','Bold');
        set(gca,'Fontsize',14);  % ,'Xlim',[0,7],'Ylim',[0,7]);
           %'XTick',[0.5,1.0,2.0,5.0,10,20],'YTick',[0.5,1.0,2.0,5.0,10,20]); % 'XScale','log','YScale','log',
        box on
        plot(gca,[0.5,4],[0.5,4],':k','Linewidth',1.5);
        str1 = {['N=',num2str(n), '; r^{2}=',num2str(mean(rsquare,'omitnan'),3),'; {\itp} = ',num2str(mean(prob,'omitnan'),3)]};
        str2 = {['Substitution Error = ',num2str(sqrt(mean(trainError,'omitnan')),3), '; Classification Error =', num2str(sqrt(mean(cvtrainError)),3)]};
        text(0.7,3.85,str1);
        text(0.7,3.65,str2);
        
        if plt_save
            print(gcf,[OutDir,'RF_regression_non-stratified_partition_Secchi_',matchup_window,'.tif'],'-dtiff','-r600');
        end
        
        figure(2);
        clf
        hp2 = bar({'rhow\_440','rhow\_480','rhow\_560','rhow\_655','rhow\_865'},mean(pred_imp,2),0.6);
        hold on
        xlabel('Predictors','Fontsize',10,'Fontname','Arial','Fontweight','Bold')
        ylabel('Predictor Importance', 'Fontsize',10,'Fontname','Arial','Fontweight','Bold');
        set(gca,'Fontsize',14);
        box on
        
        if plt_save
            print(gcf,[OutDir,'RF_regression_predictor_importance_Secchi_',matchup_window,'.tif'],'-dtiff','-r600');
        end

    case 'DO'
        if mdl_flag
            % Save RF model output
            DOMdl = Mdl;
            save([InDir,'Matchup files\RF\','RF_DO_36h_Mdl.mat'],"DOMdl");
        end

        %% Plot DO results
        
        figure(1);
        clf
        hp1 = plot(insituY,preY,'ko');
        hold on
        xlabel('RF Predicted Dissolved Oxygen (mg/L)','Fontsize',10,'Fontname','Arial','Fontweight','Bold')
        ylabel('In-Situ Dissolved Oxygen (mg/L)', 'Fontsize',10,'Fontname','Arial','Fontweight','Bold');
        set(gca,'Fontsize',14,'Xlim',[4,9],'Ylim',[4,9]);
           %'XTick',[0.5,1.0,2.0,5.0,10,20],'YTick',[0.5,1.0,2.0,5.0,10,20]); % 'XScale','log','YScale','log',
        box on
        plot(gca,[4,9],[4,9],':k','Linewidth',1.5);
        str1 = {['N=',num2str(n), '; r^{2}=',num2str(mean(rsquare,'omitnan'),3),'; {\itp} = ',num2str(mean(prob,'omitnan'),3)]};
        str2 = {['Substitution Error = ',num2str(sqrt(mean(trainError,'omitnan')),3), '; Classification Error =', num2str(sqrt(mean(cvtrainError)),3)]};
        text(4.3,8.75,str1);
        text(4.3,8.45,str2);
        
        if plt_save
            print(gcf,[OutDir,'RF_regression_non-stratified_partition_dissolved_oxygen_',matchup_window,'.tif'],'-dtiff','-r600');
        end
        
        figure(2);
        clf
        hp2 = bar({'rhow\_440','rhow\_480','rhow\_560','rhow\_655','rhow\_865'},mean(pred_imp,2),0.6);
        hold on
        xlabel('Predictors','Fontsize',10,'Fontname','Arial','Fontweight','Bold')
        ylabel('Predictor Importance', 'Fontsize',10,'Fontname','Arial','Fontweight','Bold');
        set(gca,'Fontsize',14);
        box on
        
        if plt_save
            print(gcf,[OutDir,'RF_regression_predictor_importance_dissolved_oxygen_',matchup_window,'.tif'],'-dtiff','-r600');
        end
end

disp('Completed');


