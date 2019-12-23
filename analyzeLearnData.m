clearvars;
addpath('C:\Users\yl2268\Documents\et_lotto\ASD_Pupil_analysis_scripts-master')
addpath('C:\Users\yl2268\Documents\pupil-size-master\code\dataModels')
addpath('C:\Users\yl2268\Documents\pupil-size-master\code\helperFunctions')
addpath('C:\Users\yl2268\Documents\et_lotto\RevLearn_AnalysisScript')
cd('C:\Users\yl2268\Documents\et_lotto\ProcessedData2')
%% Data processing: Step 1, screen data quality
subjs = [10:56,58:96,98:111,113:120];
sumNA = zeros(70,6);
subjNA = [];
trialNA = [];

for s = 1:length(subjs)
    load(['ET_',num2str(subjs(s))]);
    eyeData = Data.Learning.eye;
    eyeMat = eyeData{:,:};
    propNA = [];
    
    trialThr = 12;
    
    actTrial(:,s) = find(Data.Learning.rectValue == 6);
    Ori = Data.Learning.rectOri(Data.Learning.rectValue == 6);
    A = find(Data.Learning.rectOri == Ori(1));
    B = find(Data.Learning.rectOri == Ori(14));
    AacqID = setdiff(A(A<actTrial(8,s)),actTrial(:,s)); % AcquisitionA - CS+, not actualized trials
    BacqID = B(B<actTrial(8,s)); % AcquisitionB - CS-, all trials before reversal
    ArevID = A(A>actTrial(8,s));
    BrevID = setdiff(B(B>actTrial(8,s)),actTrial(:,s));
    
    AacqNum(s,1) = length(AacqID);
    BacqNum(s,1) = length(BacqID);
    ArevNum(s,1) = length(ArevID);
    BrevNum(s,1) = length(BrevID);
   
    Aacq = Data.Learning.rating(AacqID);
    Arev = Data.Learning.rating(ArevID);
    Bacq = Data.Learning.rating(BacqID);
    Brev = Data.Learning.rating(BrevID);
    
    Aacq(Aacq==0) = NaN;
    Arev(Arev==0) = NaN;
    Bacq(Bacq==0) = NaN;
    Brev(Brev==0) = NaN;    
    
    CSpAll(:,s) = [Aacq;Arev];
    CSmAll(:,s) = [Bacq;Brev];
    
    % by reversal trial number threshold
    if length(ArevID) > trialThr && length(BrevID) > trialThr
        subjLong(s,1) = subjs(s);
    AacqRa(:,s) = Data.Learning.rating(AacqID(1:14));
    BacqRa(:,s) = Data.Learning.rating(BacqID(1:14));
    ArevRa(:,s) = Data.Learning.rating(ArevID(1:trialThr));
    BrevRa(:,s) = Data.Learning.rating(BrevID(1:trialThr));
    
    AacqRt(:,s) = Data.Learning.RT(AacqID(1:14));
    BacqRt(:,s) = Data.Learning.RT(BacqID(1:14));
    ArevRt(:,s) = Data.Learning.RT(ArevID(1:trialThr));
    BrevRt(:,s) = Data.Learning.RT(BrevID(1:trialThr));
    else
    subjLong(s,1) = nan;
    AacqRa(:,s) = nan(1,14);
    BacqRa(:,s) = nan(1,14);
    ArevRa(:,s) = nan(1,trialThr);
    BrevRa(:,s) = nan(1,trialThr);
    
    AacqRt(:,s)= nan(1,14);
    BacqRt(:,s) = nan(1,14);
    ArevRt(:,s) = nan(1,trialThr);
    BrevRt(:,s)= nan(1,trialThr);
    end
    % create data model for preprocessing pipeline
    diameterUnit = 'mm';
    diameter.t_ms = [0:(1000/60):(length(eyeMat(:,1))-2)*(1000/60)]';
    diameter.L = eyeMat(2:end,5);
    diameter.R = eyeMat(2:end,8);
    zeroTime_ms = 0;
    eyeMat(:,9)=[0;diameter.t_ms];
    
    for n = 1:70
        trialData = eyeMat(eyeMat(:,1)==n,:);
        trialNA(n,s)= length(find(isnan(trialData(:,5))+isnan(trialData(:,8))>1))/length(trialData(:,1));
        % create segmentsTable
        segments.segmentStart(n,:)= trialData(1,9)/1000;
        segments.segmentEnd(n,:)= trialData(end,9)/1000;
        segments.segmentName(n,:)= {['Trial_',sprintf('%02d',n)]};
        segments.segmentSource(n,:)={['subj_',num2str(subjs(s))]};
        for m = 1:6
        segData = trialData(trialData(:,2)==m,:);
        % calculate proportional missing
        propNA(n,m) = length(find(isnan(segData(:,5))+isnan(segData(:,8)) > 1))/length(segData(:,1));
        end
    end
    segmentsTable=struct2table(segments);
    
    sumNA = sumNA+propNA;
    subjNA(s,:) = mean(propNA,1);
    %save(['convtET_',num2str(subjs(s)),'.mat'],'diameter','diameterUnit','segmentsTable','zeroTime_ms')
end
%% Behavioral data visualization
CSp = [rmmissing(AacqRa,2)',rmmissing(ArevRa,2)'];
CSm = [rmmissing(BacqRa,2)',rmmissing(BrevRa,2)'];

CSp(CSp == 0) = NaN;
CSm(CSm == 0) = NaN;

figure(1)
plot(mean(CSp,1))
errorbar([1:1:14+trialThr].',mean(CSp,1,'omitnan'), std(CSp,'omitnan')/sqrt(length(CSp)))
hold on
plot(mean(CSm,1))
errorbar([1:1:14+trialThr].',mean(CSm,1,'omitnan'), std(CSm,'omitnan')/sqrt(length(CSm)))
ylabel('Mean Ratings')
xlabel('Number of exposures')
hold on
xlim([0 15+trialThr])
ax = gca;
ax.XTick = 1:14+trialThr;
hold on
title (['Restricted: Acq 1-14, Rev 15-' num2str(14+trialThr) ' (N=' num2str(length(CSp)) ')'])
legend({'CS+', 'CS-'})
hold off

figure(2)
plot(mean(CSpAll,2))
errorbar([1:1:28].',mean(CSpAll,2,'omitnan'), std(CSpAll,0,2,'omitnan')/sqrt(length(CSpAll)))
hold on
plot(mean(CSmAll,2))
errorbar([1:1:28].',mean(CSmAll,2,'omitnan'), std(CSmAll,0,2,'omitnan')/sqrt(length(CSmAll)))
ylabel('Mean Ratings')
xlabel('Number of exposures')
hold on
xlim([0 29])
ax = gca;
ax.XTick = 1:28;
title (['All subjects (N=' num2str(length(CSpAll)) ')'])
legend({'CS+', 'CS-'})
hold off

figure (3)
boxplot([AacqNum,BacqNum,ArevNum,BrevNum],'Labels',{'Aacq','Bacq','Arev','Brev'})
AacqStats = [mean(AacqNum),std(AacqNum),range(AacqNum),median(AacqNum)];
BacqStats = [mean(BacqNum),std(BacqNum),range(BacqNum),median(BacqNum)];
ArevStats = [mean(ArevNum),std(ArevNum),range(ArevNum),median(ArevNum)];
BrevStats = [mean(BrevNum),std(BrevNum),range(BrevNum),median(BrevNum)];
%% Data Preprocessing: Step 2
% 1. Find valid subjects for calculation
% 2. zTransform, baseline correction
clearvars -except subjs trialNA subjLong
n = 1;
m = 1;

validSubj = subjs(mean(trialNA,1) <= 0.35);
%intersect(validSubj,subjLong)

actTrial = [];

for vs = 1%:length(validSubj) %54 %
    load(['ET_',num2str(validSubj(vs))]);
    load(['convtET_',num2str(validSubj(vs))]);
    eyeData = Data.Learning.eye;
    eyeMat = eyeData{:,:};
    filter = struct;
    filter.filterType = 'sgolay';
    filter.order = 3; % order of polynomial for sgolay filter?
    filter.framelen = 21; % length of window? must be odd number
    filter.clearWin = 3; % delete the n surrounding data points of a blink
    filter.velThreshold = 2; % de-blinking relative velocity threshold
    graph = true;
    
    [PupilAll,PupilL,PupilR,PupilzAll,Pupilv]=...
    combineLeftRight(diameter.L(:),diameter.R(:),diameter.t_ms(:),filter,graph);
%     propNA = [];
%     SegL = [];
%     SegR = [];
    actTrial(:,vs) = find(Data.Learning.rectValue == 6);
    Ori = Data.Learning.rectOri(Data.Learning.rectValue == 6);
    A = find(Data.Learning.rectOri == Ori(1));
    B = find(Data.Learning.rectOri == Ori(14));
    AacqID = setdiff(A(A<actTrial(8,vs)),actTrial(:,vs)); % AcquisitionA - expecting to see reward, not actualized trials
    BacqID = B(B<actTrial(8,vs)); % AcquisitionB - not expecting to see reward, all trials before reversal

    for n= 1:70
    trialData = eyeMat((eyeMat(:,1)==n),:);
    % add time column here
    %trialData(:,9) = [21/length(trialData(:,1)):21/length(trialData(:,1)):21]';
    trialData(:,9) = [0:(1000/60):(length(trialData(:,1))-1)*(1000/60)]';
 
%     [valOut,speedFiltData,devFiltData] ...
%     = rawDataFilter(trialData(:,9),trialData(:,5))
    
%     sInitial.PupilLeft = trialData(:,5);
%     sInitial.PupilRight = trialData(:,8);
    
%     [sampfiltleft,sampfiltzleft,sampinterpleft,sampdbleft,velfiltleft,veldbleft,velleft] = ...
%         pupilPrepro(trialData(:,5),trialData(:,9),filter);
%     [sampfiltright,sampfiltzright,sampinterpright,sampdbright,velfiltright,veldbright,velright] = ...
%         pupilPrepro(trialData(:,8),trialData(:,9),filter);
%     %for m = 1:6
%         SegL(:,n) = sampinterpleft((trialData(:,2)==5),:);
%         SegR(:,n) = sampinterpright((trialData(:,2)==5),:);
%         % calculate proportional missing
%         %propNA(n,m) = length(find(isnan(segData(:,5))+isnan(segData(:,8)) > 1))/length(segData(:,1));
%     %end
%
%     if sum(isnan(trialData(:,5)+trialData(:,8)))/length(trialData(:,5))*2 <= 0.5
%     check trial data quality, if more than half missing, ommit trial

% RJ's preprocessing pipeline
    [P,PL,PR,Pz,Pv]=...
    combineLeftRight(trialData(:,5),trialData(:,8),trialData(:,9),filter,graph);

%sInitial.filtered = struct;
    Pupil{n}={P};
    PupilLeft{n}={PL};
    PupilRight{n}={PR};
    Pupilz{n}={Pz};
    vel{n}={Pv};
    segOne = Pz(trialData(:,2)==1); % Segment 1, fixation point for 8 seconds
    bsWin = 1; % take the last x seconds from segment 1 as baseline
    bsline(:,n) = segOne((length(segOne)+1-bsWin*60):length(segOne),:); % should be bsWin*60 data points here
    segTwo(:,n) = Pz(trialData(:,2)==2);
%     if  (sum(isnan(bsline(:,n)))/length(bsline(:,n))) < 0.5 && (sum(isnan(segTwo(:,n)))/length(segTwo(:,n))) < 0.5
        % check baseline data quality, if more than half missing, ommit
        pSegTwo(:,n) = Pz(trialData(:,2)==2) - mean(bsline(:,n), 'omitnan'); % normalize with x second(s) before the window        
%     else
%         pSegTwo(:,n) = nan(180,1);
%         bsline(:,n) = nan(120,1);
%     end
%     sInitial.filtered = struct;
%     sInitial.filtered.Pupil(vs,n)={zeros(size(sInitial.PupilLeft))};
%     sInitial.filtered.PupilLeft(vs,n)={zeros(size(sInitial.PupilLeft))};
%     sInitial.filtered.PupilRight(vs,n)={zeros(size(sInitial.PupilLeft))};
%     sInitial.filtered.Pupilz(vs,n)={zeros(size(sInitial.PupilLeft))};
%     sInitial.filtered.vel(vs,n)={zeros(size(sInitial.PupilLeft))};
%     [sInitial.filtered.Pupil(vs,n),sInitial.filtered.PupilLeft(vs,n),sInitial.filtered.PupilRight(vs,n),sInitial.filtered.Pupilz(vs,n),sInitial.filtered.vel(vs,n)]=...
%         combineLeftRight(trialData(:,5),trialData(:,8),trialData(:,9),filter,graph);
%     else
%         segTwo(:,n) = nan(180,1);
%         bsline(:,n) = nan(120,1);
%     end
    end
    
    sumbsline = mean(bsline,1,'omitnan')';
    
    AacqEye = segTwo(:,AacqID);
    BacqEye = segTwo(:,BacqID);
    
    pAacqEye = pSegTwo(:,AacqID);
    pBacqEye = pSegTwo(:,BacqID);

    midPointA = floor(length(AacqEye(1,:))/2);
    midPointB = floor(length(BacqEye(1,:))/2);
    
    AacqE{vs} = AacqEye(:,(1:midPointA));
    AacqL{vs} = AacqEye(:,(midPointA+1:length(AacqEye(1,:))));
    BacqE{vs} = BacqEye(:,(1:midPointB));
    BacqL{vs} = BacqEye(:,(midPointB+1:length(BacqEye(1,:))));
    
    pAacqE{vs} = rmmissing(pAacqEye(:,(1:midPointA)));
    pAacqL{vs} = rmmissing(pAacqEye(:,(midPointA+1:length(AacqEye(1,:)))));
    pBacqE{vs} = rmmissing(pBacqEye(:,(1:midPointB)));
    pBacqL{vs} = rmmissing(pBacqEye(:,(midPointB+1:length(BacqEye(1,:)))));
    
    dEarly(vs) = trapz([1:180],mean(pAacqE{vs},2,'omitnan')-mean(pBacqE{vs},2,'omitnan'));
    dLate(vs) = trapz([1:180],mean(pAacqL{vs},2,'omitnan')-mean(pBacqL{vs},2,'omitnan'));
%     AvgR = (SegL(:,Data.Learning.rectValue == 6)+SegR(:,Data.Learning.rectValue == 6))/2;
%     zAvgR(:,vs) = mean(AvgR,1,'omitnan');
%     Avg = (SegL(:,Data.Learning.rectValue == 0)+SegR(:,Data.Learning.rectValue == 0))/2;
%     zAvg(:,vs) = mean(Avg,1,'omitnan');
%     figure(vs)
%     imagesc(bsline)
end

%% Pupil dynamics visualization
% Plot individual 
figure(1)
subplot(2,2,1)
for v = 1%:length(validSubj)
stdshade((AacqE{v})',0.2,'--r')
hold on
end
stdshade(cell2mat(AacqE)',0.5,'--b')
hold on
xlabel('Shape A acquisition - Early')
subplot(2,2,2)
for v = 1%:length(validSubj)
stdshade((AacqL{v})',0.2,'--c')
hold on
end
stdshade(cell2mat(AacqL)',0.5,'--b')
xlabel('Shape A acquisition - Late')
hold on
subplot(2,2,3)
for v = 1%:length(validSubj)
stdshade((BacqE{v})',0.2,'--r')
hold on
end
stdshade(cell2mat(BacqE)',0.5,'--b')
hold on
xlabel('Shape B acquisition - Early')
hold on
subplot(2,2,4)
for v = 1%:length(validSubj)
stdshade((BacqL{v})',0.2,'--c')
hold on
end
stdshade(cell2mat(BacqL)',0.5,'--b')
hold on
xlabel('Shape B acquisition - Late')
hold off


figure(2)
subplot(2,2,1)
plot(cell2mat(pAacqE))
xlabel('Shape A acquisition - Early')
hold on
subplot(2,2,2)
plot(cell2mat(pAacqL))
xlabel('Shape A acquisition - Late')
hold on
subplot(2,2,3)
plot(cell2mat(pBacqE))
xlabel('Shape B acquisition - Early')
hold on
subplot(2,2,4)
plot(cell2mat(pBacqL))
xlabel('Shape B acquisition - Late')
hold off


figure(5)
subplot(2,2,1)
imagesc(cell2mat(pAacqE)')
xlabel('Shape A acquisition - Early')
hold on
subplot(2,2,2)
imagesc(cell2mat(pAacqL)')
xlabel('Shape A acquisition - Late')
hold on
subplot(2,2,3)
imagesc(cell2mat(pBacqE)')
xlabel('Shape B acquisition - Early')
hold on
subplot(2,2,4)
imagesc(cell2mat(pBacqL)')
xlabel('Shape B acquisition - Late')
hold off
%% plot all
% Plot all subjects


figure(3)
subplot(1,2,1)
stdshade(cell2mat(AacqE)',0.5,'--r')
hold on
stdshade(cell2mat(AacqL)',0.5,'--b')
xlabel('Shape A acquisition - Early red Late blue')
hold on
subplot(1,2,2)
stdshade(cell2mat(BacqE)',0.5,'--r')
hold on
stdshade(cell2mat(BacqL)',0.5,'--b')
xlabel('Shape B acquisition - Early red Late Blue')
hold off

figure(4)
subplot(1,2,1)
stdshade(cell2mat(AacqE)',0.5,'--r')
hold on
stdshade(cell2mat(BacqE)',0.5,'--b')
xlabel('Early Acq - CS+ red CS- blue')
hold on
subplot(1,2,2)
stdshade(cell2mat(AacqL)',0.5,'--r')
hold on
stdshade(cell2mat(BacqL)',0.5,'--b')
xlabel('Late Acq - CS+ red CS- blue')
hold off

figure(5)
subplot(1,2,1)
stdshade(cell2mat(pAacqE)',0.5,'--r')
hold on
stdshade(cell2mat(pBacqE)',0.5,'--b')
xlabel('Early Acq - CS+ red CS- blue')
ylim([-30000 30000])
hold on
subplot(1,2,2)
stdshade(cell2mat(pAacqL)',0.5,'--r')
hold on
stdshade(cell2mat(pBacqL)',0.5,'--b')
xlabel('Late Acq - CS+ red CS- blue')
ylim([-30000 30000])
hold off

%% Subject early/late contrast
for sj =20 %1:length(validSubj)
figure(sj)
subplot(1,2,1)
for v = sj%1:length(validSubj)
stdshade((pAacqE{v})',0.2,'--m')
hold on
stdshade((pBacqE{v})',0.2,'--c')
end
% stdshade(cell2mat(pAacqE)',0.7,'--r')
% hold on
% stdshade(cell2mat(pBacqE)',0.7,'--b')
hold on
xlabel('Early Acquisition')
%ylim([-1000 1000])
% hold on
subplot(1,2,2)
for v = sj%1:length(validSubj)
stdshade((pAacqL{v})',0.2,'--m')
hold on
stdshade((pBacqL{v})',0.2,'--c')
hold on
end
% stdshade(cell2mat(pAacqL)',0.5,'--r')
% hold on
% stdshade(cell2mat(pBacqL)',0.5,'--b')
hold on
xlabel('Late Acquisition')
hold on
legends({'CS+','CS-'})
%ylim([-1000 1000])
hold off
end
% figure(3)
% subplot(2,2,1)
% imagesc(cell2mat(AacqE)')
% xlabel('Shape A acquisition - Early')
% hold on
% subplot(2,2,2)
% imagesc(cell2mat(AacqL)')
% xlabel('Shape A acquisition - Late')
% hold on
% subplot(2,2,3)
% imagesc(cell2mat(BacqE)')
% xlabel('Shape B acquisition - Early')
% hold on
% subplot(2,2,4)
% imagesc(cell2mat(BacqL)')
% xlabel('Shape B acquisition - Late')
% holf off

% %% Data visualization
% avgNA = sumNA / length(subjs);
% figure(1)
% imagesc(avgNA)
% hold on
% title('Average missing, segment by trial')
% xlabel('Segment')
% ylabel('Trial Number')
% hold off
% 
% figure(2)
% boxplot(subjNA)
% hold on
% title('Proportion missing by segment')
% xlabel('Segment')
% ylabel('Trial Number')
% hold off
% 
% figure(3)
% imagesc(trialNA)
% trsTrial = trialNA(:,(mean(trialNA,1) < 0.35));
% imagesc(trsTrial)
% 
% STD = std(zscore(zAvg),0,2,'omitnan')/sqrt(length(zAvg(1,:)));
% MEAN = mean(zscore(zAvg),2,'omitnan');
% 
% figure(4)
% x = 1:1:120;
% y = MEAN';
% % I create some yu and yl here, for the example
% yu = y+STD;
% yl = y-STD;
% fill([x fliplr(x)], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none')
% hold all
% plot(x,y)
% 
% STD = std(zscore(zAvgR),0,2,'omitnan')/sqrt(length(zAvgR(1,:)));
% MEAN = mean(zscore(zAvgR),2,'omitnan');
% 
% figure(5)
% x = 1:1:120;
% y = MEAN';
% % I create some yu and yl here, for the example
% yu = y+STD;
% yl = y-STD;
% fill([x fliplr(x)], [yu fliplr(yl)], [.9 .9 .9], 'linestyle', 'none')
% hold all
% plot(x,y)



