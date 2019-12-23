clearvars;
cd('C:\Users\yl2268\Documents\et_lotto\ProcessedData2')
%% for all subjects
subj = [10:96,98:120];
for s = 1:length(subj)
    % load Data from each subject .mat file within the folder
    load(['ET_',num2str(subj(s))]);
    % add proportion FOSD violation
    ratioAll(s) = Data.Lotto.fivRatio;
end
% mean and std of FOSD violation proportion for all
Sum.avgFOSDall = {mean(ratioAll),std(ratioAll)};
%% for valid subjects
validSubj = subj(ratioAll < 0.5);
n = 1;
m = 1;
for i = 1:length(validSubj)
    % load Data from each subject .mat file within the folder
    load(['ET_',num2str(validSubj(i))]);
    % add proportion FOSD violation
    ratioValid(i) = Data.Lotto.fivRatio;
    % add proportion of lottory choice (for all valid and for rep. trials)
    lottoPropValid(i) = Data.Lotto.lottoPropAll;
    lottoPropRep(i) = Data.Lotto.lottoPropRep;
    % add alpha and gamma values (for all valid and for rep. trials)
    fitValid(i,:)=Data.Lotto.fitAll;
    fitRep(i,:)=Data.Lotto.fitRep;
    % add FBI mean for each subject
    FBImeanAll(i) = Data.Lotto.FBImean;
    Type(i) = Data.Lotto.Type;
    %% based on group - All
    if Data.Lotto.Type ==1
    FBImeanOne(n)= Data.Lotto.FBImean;
    alphaOne(1,n) = Data.Lotto.fitAll(1);
%     alphaOne(2,n) = validSubj(i);
%     alphaOne(3,n) = 1;
    
    gammaOne(1,n) = Data.Lotto.fitAll(2);
%     gammaOne(2,n) = validSubj(i);
%     gammaOne(3,n) = 1;
    
    pweightOne(1,n) = Data.Lotto.fitAll(3);
% pweightOne(2,n) = validSubj(i);     
% pweightOne(3,n) = 1;
    alphaOneR(1,n) = Data.Lotto.fitRep(1);
% alphaOneR(2,n) = validSubj(i);
% alphaOneR(3,n) = 1;
    
    gammaOneR(1,n) = Data.Lotto.fitRep(2);
% gammaOneR(2,n) = validSubj(i);
% gammaOneR(3,n) = 1;
    
    pweightOneR(1,n) = Data.Lotto.fitRep(3);
% pweightOneR(2,n) = validSubj(i);
% pweightOneR(3,n) = 1;
    
    n = n+1;
    else
    FBImeanTwo(m)= Data.Lotto.FBImean;
    alphaTwo(1,m) = Data.Lotto.fitAll(1);
%     alphaTwo(2,m) = validSubj(i);
%     alphaTwo(3,m) = 2;
    
    gammaTwo(1,m) = Data.Lotto.fitAll(2);
%     gammaTwo(2,m) = validSubj(i);
%     gammaTwo(3,m) = 2;
    
    pweightTwo(1,m) = Data.Lotto.fitAll(3);
%      pweightTwo(2,m) = validSubj(i);
%      pweightTwo(3,m) = 2;
    alphaTwoR(1,m) = Data.Lotto.fitRep(1);
% alphaTwoR(2,m) = validSubj(i);
% alphaTwoR(3,m) = 2;
    
    gammaTwoR(1,m) = Data.Lotto.fitRep(2);
% gammaTwoR(2,m) = validSubj(i);
% gammaTwoR(3,m) = 2;
    
    pweightTwoR(1,m) = Data.Lotto.fitRep(3);
% pweightTwoR(2,m) = validSubj(i);
% pweightTwoR(3,m) = 2;
    
    m = m+1;
    end
    %% based on group - rep trials
%     if Data.Lotto.Type ==1
%     alphaOneR(1,n) = Data.Lotto.fitRep(1);
%     alphaOneR(2,n) = validSubj(i);
%     alphaOneR(3,n) = 1;
%     
%     gammaOneR(1,n) = Data.Lotto.fitRep(2);
%     gammaOneR(2,n) = validSubj(i);
%     gammaOneR(3,n) = 1;
%     
%     pweightOneR(1,n) = Data.Lotto.fitRep(3);
%      pweightOneR(2,n) = validSubj(i);
%      pweightOneR(3,n) = 1;
%     
%     n = n+1;
%     else
%     alphaTwoR(1,m) = Data.Lotto.fitRep(1);
% %     alphaTwoR(2,m) = validSubj(i);
% %     alphaTwoR(3,m) = 2;
%     
%     gammaTwoR(1,m) = Data.Lotto.fitRep(2);
% %     gammaTwoR(2,m) = validSubj(i);
% %     gammaTwoR(3,m) = 2;
%     
%     pweightTwoR(1,m) = Data.Lotto.fitRep(3);
%      pweightTwoR(2,m) = validSubj(i);
%      pweightTwoR(3,m) = 2;
%     
%     m = m+1;
%     end
end
fitValid(:,3) = NaN(length(fitValid),1);
fitValid(find(fitValid(:,1) > 1),3) = zeros;
fitValid(find(fitValid(:,1) < 1),3) = ones;
% mean and std of FOSD violation proportion for valid trials
Sum.avgFOSDvalid = {mean(ratioValid),std(ratioValid)}; 
% mean and std of lottory proportion for all trials
Sum.LottoValid= {mean(lottoPropValid), std(lottoPropValid)};
% mean and std of lottory proportion for rep. trials
Sum.LottoRep = {mean(lottoPropRep),std(lottoPropRep)};
% mean and std of alpha for all trials, Pic 4
Sum.fitValidAlpha = {mean(fitValid(:,1)),std(fitValid(:,1))};
% mean and std of alpha proportion for rep. trials
Sum.fitRepAlpha = {mean(fitRep(:,1)),std(fitRep(:,1))};
% mean and std of gamma for all trials
Sum.fitValidGamma = {mean(fitValid(:,2)),std(fitValid(:,2))};
% mean and std of gamma proportion for rep. trials
Sum.fitRepGamma = {mean(fitRep(:,2)),std(fitRep(:,2))};
% mean and std of probability weight for rep. trials
Sum.fitReppweight = {mean(fitRep(:,3)),std(fitRep(:,3))};
% mean and std of FBI for all
Sum.FBImean = {mean(FBImeanAll),std(FBImeanAll)};

Sum.Group1FBI = {mean(FBImeanOne),std(FBImeanOne)};
Sum.Group1Alpha = {mean(alphaOne),std(alphaOne)};
Sum.Group1Gamma = {mean(gammaOne),std(gammaOne)};
Sum.Group1pweight = {mean(pweightOne),std(pweightOne)};
Sum.Group2FBI = {mean(FBImeanTwo),std(FBImeanTwo)};
Sum.Group2Alpha = {mean(alphaTwo),std(alphaTwo)};
Sum.Group2Gamma = {mean(gammaTwo),std(gammaTwo)};
Sum.Group2pweight = {mean(pweightTwo),std(pweightTwo)};
SumFitAll = struct2table(Sum);

AlphaAll = [alphaOne';alphaTwo'];
GammaAll = [gammaOne';gammaTwo'];
pweightAll = [pweightOne';pweightTwo'];

AlphaOvl = [alphaOneR';alphaTwoR'];
GammaOvl = [gammaOneR';gammaTwoR'];
pweightOvl = [pweightOneR';pweightTwoR'];

%% Other tests
% Spearman correlation
[RHO,PVAL] = corr([FBImeanOne(:);FBImeanTwo(:)],[pweightOne(:);pweightTwo(:)],'Type','Spearman');
% FBI
% Wilcoxon rank sum test
WRST.FBI = ranksum(FBImeanOne, FBImeanTwo);
% Randomiaztion test
RmedDiff = median(FBImeanOne) - median(FBImeanTwo);
[FBIperm,FBI25,FBI75,FBIpVal] = permutation(10000,FBImeanOne,FBImeanTwo,RmedDiff);
% sum(RmedDiff<=abs(FBIperm))
% medProp = sum(FBIperm >= RmedDiff)/10000;
% % plot distribution
% figure
% hist(FBIperm,100)

% Risk attitute
WRST.Risk = ranksum(alphaOne,alphaTwo);
RmedDiffRisk = median(alphaOne) - median(alphaTwo);
[Riskperm,Risk25,Risk75,RiskpVal] = permutation(10000,alphaOne,alphaTwo,RmedDiffRisk);
% figure
% hist(Riskperm,100)

% Probability Weighting
WRST.pweight = ranksum(pweightOne,pweightTwo);
RmedDiffpweight = median(pweightOne) - median(pweightTwo);
[pweightperm,pweight25,pweight75,pweightpVal] = permutation(10000,pweightOne,pweightTwo,RmedDiffpweight);

% Repeated trials Probability Weighting
WRST.pweightR = ranksum(pweightOneR(1,:),pweightTwoR(1,:));
RmedDiffpweightR = median(pweightOneR(1,:)) - median(pweightTwoR(1,:));
[pweightpermR,pweight25R,pweight75R,pweightpValR] = permutation(10000,pweightOneR(1,:),pweightTwoR(1,:),RmedDiffpweightR);

% confidence interval

%%
% %% Plot
% % Pic 1. mean and sem of FBI by group
% group = categorical({'Val>Prob','Prob>Val'});
% meanFBI = [Sum.Group1FBI{1}, Sum.Group2FBI{1}];
% semFBI = [Sum.Group1FBI{2}/sqrt(length(FBImeanOne)),Sum.Group2FBI{2}/sqrt(length(FBImeanTwo))];
% subplot(2,2,1)
% bar(group,meanFBI)
% xlabel('Subject Group')
% ylabel('FBI')
% title('Plot 1: Mean FBI by group')
% hold on 
% errorbar(group,meanFBI,semFBI,'k','linestyle','none')
% hold off
% 
% % Pic 2. mean and sem of lottery prop for overlapping choice set
% meanLottoPropRep = [Sum.LottoRep{1}];
% semLottoPropRep = [Sum.LottoRep{2}/sqrt(length(lottoPropRep))];
% subplot(2,2,2)
% bar(categorical({'All Subj'}),meanLottoPropRep)
% xlabel('Subject Group')
% ylabel('Proportion of choosing lottery (OVLP)')
% title('Plot 2: mean prop. chose lottery of OVLP set')
% hold on
% errorbar(categorical({'All Subj'}),meanLottoPropRep,semLottoPropRep,'k','linestyle','none')
% hold off
% 
% % Pic 3. mean and sem of alpha by group
% meanAlpha = [Sum.Group1Alpha{1}, Sum.Group2Alpha{1}];
% semAlpha = [Sum.Group1Alpha{2}/sqrt(length(alphaOne)),Sum.Group2Alpha{2}/sqrt(length(alphaTwo))];
% subplot(2,2,3)
% bar(group,meanAlpha)
% xlabel('Subject Group')
% ylabel('Alpha')
% title('Plot 3: mean alpha by group')
% hold on
% errorbar(group,meanAlpha,semAlpha,'k','linestyle','none')
% hold off
% 
% % Pic 4. mean and sem of alpha for full set
% subplot(2,2,4)
% bar(categorical({'All Subj'}),Sum.fitValidAlpha{1})
% xlabel('Subject Group')
% ylabel('Alpha')
% title('Plot 4: mean alpha for all subjects')
% hold on
% errorbar(categorical({'All Subj'}),Sum.fitValidAlpha{1},Sum.fitValidAlpha{2}/sqrt(length(fitValid(:,1))),'k','linestyle','none')
% hold off
% 
% %% Save pic
% % savefig(['summary_subj',num2str(min(subj)),'to',num2str(max(subj)),'.fig'])
