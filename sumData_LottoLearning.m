clearvars;
cd('C:\Users\yl2268\Documents\et_lotto\ProcessedData2')
%% for all subjects
subj = [10:96,98:120];
for s = 1:length(subj)
    % load Data from each subject .mat file within the folder
    %print(subj(s))
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
nn = 1;
mm = 1;
for i = 1:length(validSubj)
    %% Lotto
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
    %% based on group
    if Data.Lotto.Type ==1
    FBImeanOne(n)= Data.Lotto.FBImean;
    alphaOne(n) = Data.Lotto.fitAll(1);
    n = n+1;
    else
    FBImeanTwo(m)= Data.Lotto.FBImean;
    alphaTwo(m) = Data.Lotto.fitAll(1);
    m = m+1;
    end
    
    %% Learning
    % load Data from each subject .mat file within the folder
    % load(['ET_',num2str(validSubj(i))]);
    % identify CS order
    Ori0Ind = intersect(find(Data.Learning.rectValue == 6),find(Data.Learning.rectOri == 0));
    Ori45Ind = intersect(find(Data.Learning.rectValue == 6), find(Data.Learning.rectOri == 45));
    %if Data.Learning.rectValue(1) ~=6
    if max(Ori0Ind) <= 35
        Group(i) = 0; % reward is associated with + at start off
        plus.rating(:,nn) = Data.Learning.rating; % add rating history of Group +
        plus.reward(:,nn) = Data.Learning.rectValue; % add reward history of Group +
        NonRefIndPlus = find(Data.Learning.rectOri == 0);
        NonRefIndX = find(Data.Learning.rectOri == 45);
        P.acqPlusCS1(:,nn) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndPlus(NonRefIndPlus<=35))); % 1
        P.revPlusCS0(:,nn) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndPlus(NonRefIndPlus>35))); % 2
        P.revXCS1(:,nn) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndX(NonRefIndX>35))); % 4
        P.acqXCS0(:,nn) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndX(NonRefIndX<=35))); % 3
        P.LottovLearn(nn) = Data.Lotto.lottoPropRep;
        nn = nn + 1;
    else
        Group(i) = 1; % reward is associated with x at start off
        x.rating(:,mm) = Data.Learning.rating; % add rating history of Group x 
        x.reward(:,mm) = Data.Learning.rectValue; % add rating history of Group x 
        NonRefIndX = find(Data.Learning.rectOri == 45);
        NonRefIndPlus = find(Data.Learning.rectOri == 0);
        X.revPlusCS1(:,mm) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndPlus(NonRefIndPlus>35))); % 8
        X.acqPlusCS0(:,mm) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndPlus(NonRefIndPlus<=35))); % 7
        X.acqXCS1(:,mm) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndX(NonRefIndX<=35))); % 5
        X.revXCS0(:,mm) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndX(NonRefIndX>35))); % 6
        X.LottovLearn(mm) = Data.Lotto.lottoPropRep;
        mm = mm + 1; 
    end
    %end
end


%% Lotto summary
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
% mean and std of FBI for all
Sum.FBImean = {mean(FBImeanAll),std(FBImeanAll)};

Sum.Group1FBI = {mean(FBImeanOne),std(FBImeanOne)};
Sum.Group1Alpha = {mean(alphaOne),std(alphaOne)};
Sum.Group2FBI = {mean(FBImeanTwo),std(FBImeanTwo)};
Sum.Group2Alpha = {mean(alphaTwo),std(alphaTwo)};

SumFitAll = struct2table(Sum);
%% Lotto Plot
% Pic 1. mean and sem of FBI by group
group = categorical({'Val>Prob','Prob>Val'});
meanFBI = [Sum.Group1FBI{1}, Sum.Group2FBI{1}];
semFBI = [Sum.Group1FBI{2}/sqrt(length(FBImeanOne)),Sum.Group2FBI{2}/sqrt(length(FBImeanTwo))];
subplot(2,2,1)
bar(group,meanFBI)
xlabel('Subject Group')
ylabel('FBI')
title('Plot 1: Mean FBI by group')
hold on 
errorbar(group,meanFBI,semFBI,'k','linestyle','none')
hold off

% Pic 2. mean and sem of lottery prop for overlapping choice set
meanLottoPropRep = [Sum.LottoRep{1}];
semLottoPropRep = [Sum.LottoRep{2}/sqrt(length(lottoPropRep))];
subplot(2,2,2)
bar(categorical({'All Subj'}),meanLottoPropRep)
xlabel('Subject Group')
ylabel('Proportion of choosing lottery (OVLP)')
title('Plot 2: mean prop. chose lottery of OVLP set')
hold on
errorbar(categorical({'All Subj'}),meanLottoPropRep,semLottoPropRep,'k','linestyle','none')
hold off

% Pic 3. mean and sem of alpha by group
meanAlpha = [Sum.Group1Alpha{1}, Sum.Group2Alpha{1}];
semAlpha = [Sum.Group1Alpha{2}/sqrt(length(alphaOne)),Sum.Group2Alpha{2}/sqrt(length(alphaTwo))];
subplot(2,2,3)
bar(group,meanAlpha)
xlabel('Subject Group')
ylabel('Alpha')
title('Plot 3: mean alpha by group')
hold on
errorbar(group,meanAlpha,semAlpha,'k','linestyle','none')
hold off

% Pic 4. mean and sem of alpha for full set
subplot(2,2,4)
bar(categorical({'All Subj'}),Sum.fitValidAlpha{1})
xlabel('Subject Group')
ylabel('Alpha')
title('Plot 4: mean alpha for all subjects')
hold on
errorbar(categorical({'All Subj'}),Sum.fitValidAlpha{1},Sum.fitValidAlpha{2}/sqrt(length(fitValid(:,1))),'k','linestyle','none')
hold off

%% Lotto Save pic
savefig(['summary_subj',num2str(min(subj)),'to',num2str(max(subj)),'.fig'])

%% Learning
% Replace zeros (missing rating with NaN)
% P.rating(plus.rating == 0) = NaN;
P.acqPlusCS1(P.acqPlusCS1 == 0) = NaN;
P.revPlusCS0(P.revPlusCS0 == 0) = NaN;
P.revXCS1(P.revXCS1 == 0) = NaN;
P.acqXCS0(P.acqXCS0 == 0) = NaN;
%plus.summary = struct2table(P);

% X.rating(x.rating == 0) = NaN;
X.acqXCS1(X.acqXCS1 == 0) = NaN;
X.revXCS0(X.revXCS0 == 0) = NaN;
X.acqPlusCS0(X.acqPlusCS0 == 0) = NaN;
X.revPlusCS1(X.revPlusCS1 == 0) = NaN;
%x.summary = struct2table(X);

% Concatenate by group
Aacq = [P.acqPlusCS1,X.acqXCS1]; % 1 + 5
Arev = [P.revPlusCS0,X.revXCS0]; % 2 + 6
Bacq = [P.acqXCS0,X.acqPlusCS0]; % 3 + 7
Brev = [P.revXCS1,X.revPlusCS1]; % 4 + 8

% Find mean and std by shape A and B
Aavg = nanmean([Aacq;Arev],2);
Astd = nanstd([Aacq;Arev],0,2);
Asem = Astd/sqrt(length(Aacq));

Bavg = nanmean([Bacq;Brev],2);
Bstd = nanstd([Bacq;Brev],0,2);
Bsem = Bstd/sqrt(length(Bacq));

% Divide by eary(the first half) and late(the second half)
% Acquisition
AcqE.A = [mean(nanmean(Aacq(1:7,:),1)),std(nanmean(Aacq(1:7,:),1))];
AcqE.A(3) = AcqE.A(2)/sqrt(length(validSubj));
AcqL.A = [mean(nanmean(Aacq(8:14,:),1)),std(nanmean(Aacq(8:14,:),1))];
AcqL.A(3) = AcqL.A(2)/sqrt(length(validSubj));

AcqE.B = [mean(nanmean(Bacq(1:7,:),1)),std(nanmean(Bacq(1:7,:),1))];
AcqE.B(3) = AcqE.B(2)/sqrt(length(validSubj));
AcqL.B = [mean(nanmean(Bacq(8:14,:),1)),std(nanmean(Bacq(8:14,:),1))];
AcqL.B(3) = AcqL.B(2)/sqrt(length(validSubj));

% Reversal
RevE.A = [mean(nanmean(Arev(1:7,:),1)),std(nanmean(Arev(1:7,:),1))];
RevE.A(3) = RevE.A(2)/sqrt(length(validSubj));
RevL.A = [mean(nanmean(Arev(8:14,:),1)),std(nanmean(Arev(8:14,:),1))];
RevL.A(3) = RevL.A(2)/sqrt(length(validSubj));

RevE.B = [mean(nanmean(Brev(1:7,:),1)),std(nanmean(Brev(1:7,:),1))];
RevE.B(3) = RevE.B(2)/sqrt(length(validSubj));
RevL.B = [mean(nanmean(Brev(8:14,:),1)),std(nanmean(Brev(8:14,:),1))];
RevL.B(3) = RevL.B(2)/sqrt(length(validSubj));

%% Summarize data - Lotto v Learning
% Acq Late
LvLAcqL(:,1) = lottoPropRep;
LvLAcqL(:,2) = nanmean(Aacq(8:14,:),1);
LvLAcqL(:,3) = nanmean(Bacq(8:14,:),1);

% Rev Late
LvLRevL(:,1) = lottoPropRep;
LvLRevL(:,2) = nanmean(Arev(8:14,:),1);
LvLRevL(:,3) = nanmean(Brev(8:14,:),1);

%% Plot - difference in rating of A and B
x = [P.LottovLearn,X.LottovLearn];
y1 = nanmean(Aacq(8:14,:),1) - nanmean(Bacq(8:14,:),1);
subplot(2,1,1)
scatter(x,y1)
hold on
legend({'Shape A - Shape B'})
ylabel('Difference of Mean Ratings')
xlabel('Probability of choosing lottery')
title('Late Acquisition')

% Reversal
y2 = nanmean(Arev(8:14,:),1) - nanmean(Brev(8:14,:),1);
subplot(2,1,2)
scatter(x,y2)
hold on
legend({'Shape A - Shape B'})
ylabel('Difference of Mean Ratings')
xlabel('Probability of choosing lottery')
title('Late Reversal')

%% Plot - Lotto v Learning
% Acquisition
b1 = LvLAcqL(:,1)\LvLAcqL(:,2);
yCalc1 = b1*LvLAcqL(:,1);
b2 = LvLAcqL(:,1)\LvLAcqL(:,3);
yCalc2 = b2*LvLAcqL(:,1);
subplot(2,1,1)
scatter(LvLAcqL(:,1),LvLAcqL(:,2))
hold on
scatter(LvLAcqL(:,1),LvLAcqL(:,3))
legend({'Shape A','Shape B'})
hold on
SlopeA = plot(LvLAcqL(:,1),yCalc1)
SlopeB = plot(LvLAcqL(:,1),yCalc2,'--')
ylabel('Mean Ratings at Late Acquisition')
xlabel('Probability of choosing lottery')

% Reversal
b3 = LvLRevL(:,1)\LvLRevL(:,2);
yCalc3 = b3*LvLRevL(:,1);
b4 = LvLRevL(:,1)\LvLRevL(:,3);
yCalc4 = b4*LvLRevL(:,1);

subplot(2,1,2)
scatter(LvLRevL(:,1),LvLRevL(:,2))
hold on
scatter(LvLRevL(:,1),LvLRevL(:,3))
legend({'Shape A','Shape B'})
hold on
SlopeA = plot(LvLRevL(:,1),yCalc3)
SlopeB = plot(LvLRevL(:,1),yCalc4,'--')
ylabel('Mean Ratings at Late Reversal')
xlabel('Probability of choosing lottery')