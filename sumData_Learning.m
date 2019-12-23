clearvars;
cd('C:\Users\yl2268\Documents\et_lotto\ProcessedData')
%% for all subjects
subj = 10:45;
n = 1;
m = 1;
for s = 1:length(subj)
    % load Data from each subject .mat file within the folder
    load(['ET_',num2str(subj(s))]);
    % identify CS order
    Ori0Ind = intersect(find(Data.Learning.rectValue == 6),find(Data.Learning.rectOri == 0));
    Ori45Ind = intersect(find(Data.Learning.rectValue == 6), find(Data.Learning.rectOri == 45));
    if Data.Learning.rectValue(1) ~=6
    if max(Ori0Ind) <= 35
        Group(s) = 0; % reward is associated with + at start off
        plus.rating(:,n) = Data.Learning.rating; % add rating history of Group +
        plus.reward(:,n) = Data.Learning.rectValue; % add reward history of Group +
        NonRefIndPlus = find(Data.Learning.rectOri == 0);
        NonRefIndX = find(Data.Learning.rectOri == 45);
        P.acqPlusCS1(:,n) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndPlus(NonRefIndPlus<=35))); % 1
        P.revPlusCS0(:,n) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndPlus(NonRefIndPlus>35))); % 2
        P.revXCS1(:,n) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndX(NonRefIndX>35))); % 4
        P.acqXCS0(:,n) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndX(NonRefIndX<=35))); % 3
        n = n + 1;
    else
        Group(s) = 1; % reward is associated with x at start off
        x.rating(:,m) = Data.Learning.rating; % add rating history of Group x 
        x.reward(:,m) = Data.Learning.rectValue; % add rating history of Group x 
        NonRefIndX = find(Data.Learning.rectOri == 45);
        NonRefIndPlus = find(Data.Learning.rectOri == 0);
        X.revPlusCS1(:,m) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndPlus(NonRefIndPlus>35))); % 8
        X.acqPlusCS0(:,m) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndPlus(NonRefIndPlus<=35))); % 7
        X.acqXCS1(:,m) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndX(NonRefIndX<=35))); % 5
        X.revXCS0(:,m) = Data.Learning.rating(intersect(find(Data.Learning.rectValue ~= 6),NonRefIndX(NonRefIndX>35))); % 6
        m = m + 1; 
    end
    end
end
% Replace zeros (missing rating with NaN)
% P.rating(plus.rating == 0) = NaN;
P.acqPlusCS1(P.acqPlusCS1 == 0) = NaN;
P.revPlusCS0(P.revPlusCS0 == 0) = NaN;
P.revXCS1(P.revXCS1 == 0) = NaN;
P.acqXCS0(P.acqXCS0 == 0) = NaN;
plus.summary = struct2table(P);

% X.rating(x.rating == 0) = NaN;
X.acqXCS1(X.acqXCS1 == 0) = NaN;
X.revXCS0(X.revXCS0 == 0) = NaN;
X.acqPlusCS0(X.acqPlusCS0 == 0) = NaN;
X.revPlusCS1(X.revPlusCS1 == 0) = NaN;
x.summary = struct2table(X);

% Concatenate by group
Aacq = [P.acqPlusCS1,X.acqXCS1]; % 1 + 5
Arev = [P.revPlusCS0,X.revXCS0]; % 2 + 6
Bacq = [P.acqXCS0,X.acqPlusCS0]; % 3 + 7
Brev = [P.revXCS1,X.revPlusCS1]; % 4 + 8
% test
[h1,p1]= ttest(nanmean(Aacq(1:7,:)),nanmean(Bacq(1:7,:)));
[h2,p2]= ttest(nanmean(Aacq(8:14,:)),nanmean(Bacq(8:14,:)));
[h3,p3]= ttest(nanmean(Arev(1:7,:)),nanmean(Brev(1:7,:)));
[h4,p4]= ttest(nanmean(Arev(8:14,:)),nanmean(Brev(8:14,:)));
pValue = num2str([p1;p2;p3;p4]);
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
AcqE.A(3) = AcqE.A(2)/sqrt(s);
AcqL.A = [mean(nanmean(Aacq(8:14,:),1)),std(nanmean(Aacq(8:14,:),1))];
AcqL.A(3) = AcqL.A(2)/sqrt(s);

AcqE.B = [mean(nanmean(Bacq(1:7,:),1)),std(nanmean(Bacq(1:7,:),1))];
AcqE.B(3) = AcqE.B(2)/sqrt(s);
AcqL.B = [mean(nanmean(Bacq(8:14,:),1)),std(nanmean(Bacq(8:14,:),1))];
AcqL.B(3) = AcqL.B(2)/sqrt(s);

% Reversal
RevE.A = [mean(nanmean(Arev(1:7,:),1)),std(nanmean(Arev(1:7,:),1))];
RevE.A(3) = RevE.A(2)/sqrt(s);
RevL.A = [mean(nanmean(Arev(8:14,:),1)),std(nanmean(Arev(8:14,:),1))];
RevL.A(3) = RevL.A(2)/sqrt(s);

RevE.B = [mean(nanmean(Brev(1:7,:),1)),std(nanmean(Brev(1:7,:),1))];
RevE.B(3) = RevE.B(2)/sqrt(s);
RevL.B = [mean(nanmean(Brev(8:14,:),1)),std(nanmean(Brev(8:14,:),1))];
RevL.B(3) = RevL.B(2)/sqrt(s);

%% Plot
% Plot mean ratings by number of exposures
subplot(2,1,1)
plot([1:28].',Aavg,[1:28].',Bavg)
errorbar([1:28].',Aavg, Asem)
hold on
errorbar([1:28].',Bavg, Bsem)
lgd = legend({'Shape A','Shape B'})
ylabel('Mean Ratings')
xlabel('Number of exposures')
xlim([0 30])
ax = gca;
ax.XTick = 1:28;

% Plot bar graph of mean ratings by groups
Avg = [AcqE.A(1),AcqE.B(1); AcqL.A(1),AcqL.B(1); RevE.A(1),RevE.B(1); RevL.A(1),RevL.B(1)];
Sem = [AcqE.A(3),AcqE.B(3); AcqL.A(3),AcqL.B(3); RevE.A(3),RevE.B(3); RevL.A(3),RevL.B(3)];
group = {'Early Acq','Late Acq','Early Rev','Late Rev'};
subplot(2,1,2)
hBar = bar(Avg,1);
hold on 
xBar = cell2mat(get(hBar,'XData')).' + [hBar.XOffset];
hEB = errorbar(xBar,Avg,Sem,'k.');
set(gca,'xtick',[1:length(group)],'xticklabel',group)
ylabel('Mean Ratings')
legend({'Shape A', 'Shape B'})