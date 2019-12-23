function processData(subjs)
root = 'C:\Users\yl2268\Documents\et_lotto\RawData\';
savedir = 'C:\Users\yl2268\Documents\et_lotto\ProcessedData2';
for s = 1:length(subjs)
    %% Find the lotto experiment .csv of this subject
    subj = num2str(subjs(s));
    Data.filename = fullfile(savedir,['ET_' num2str(subj)]);
    LottoCSV = dir([root,'**\ETLotto*_',subj,'.csv']);
    for n = 1:length(LottoCSV)
        name = LottoCSV(n).name;
        idx = strfind(name,'Practice');
        if isempty(idx)
            Lotto = LottoCSV(n);
        end
    end
     %% Find the last learning .csv if multiple
    %Learning folder might have more than one csv data files
    LearningCSV = dir([root,'**\ETLearning*',subj,'.csv']);
    if length(LearningCSV) > 1
        [val,idx] = max([LearningCSV.bytes]);
        Learning = LearningCSV(idx);
    else
        Learning = LearningCSV;
    end
    %% Convert table to structure 
    Data.Lotto = table2struct(readtable([Lotto.folder,'\',Lotto.name]),'ToScalar',true);
    Data.Learning = table2struct(readtable([Learning.folder,'\',Learning.name]),'ToScalar',true);
    %% Add eyetracking data
    etNames = {'TrialNum','SegNum','gazeXR','gazeYR','pupilR','gazeXL','gazeYL','pupilL'};
    % Lotto
    LottoDir = dir([root,'**\ETLotto_',name(9:19),'*.mat']); % find the .mat file for this subj
    LottoLoad = load([LottoDir.folder,'\',LottoDir.name]); % load file
    Data.Lotto.eye = array2table(LottoLoad.etData(:,[1,2,14,15,24,33,34,43]),'VariableNames',...
    etNames);
    % Learning
    LearningDir = dir([root,'**\ETLearning_',Learning.name(12:22),'*.mat']); % find the .mat file for this subj
    LearningLoad = load([LearningDir.folder,'\',LearningDir.name]); % load file
    Data.Learning.eye = array2table(LearningLoad.etData(:,[1,2,14,15,24,33,34,43]),'VariableNames',...
        etNames);
    %% Calculate proportion of choosing lottery when val = $5
    fivInd = find(Data.Lotto.lottoVal == 5);
    Data.Lotto.fivRatio = sum(Data.Lotto.choseLotto(fivInd))/length(fivInd);
    %% Model based & model free calculation - for all trials
    choseLottery = Data.Lotto.choseLotto;
    refVal = ones(size(choseLottery))*5;
    lottoVal = Data.Lotto.lottoVal;
    refProb = ones(size(choseLottery));
    lottoProb = Data.Lotto.lottoProb/100;
    b0 = [1,1,1]; % initial estimate of alpha and gamma
    % choseLotto,refVal,lottoVal,refProb,lottoProb,b0
    fittingResults = fitData_SU(choseLottery,refVal,lottoVal,refProb,lottoProb,b0);
    alpha=fittingResults.b(1);
    gamma=fittingResults.b(2);
    pweight=fittingResults.b(3);
    Data.Lotto.fitAll=[alpha,gamma,pweight];
    % proportion of lottery choices
    Data.Lotto.lottoPropAll = sum(Data.Lotto.choseLotto)/length(Data.Lotto.choseLotto);
    %% Model based & model free calculation - for overlapped trials
    % repeat model-based and model-free for overlapping choice set
    if length(unique(Data.Lotto.lottoVal)) > length(unique(Data.Lotto.lottoProb))
        Data.Lotto.Type = 1; % group 1 = no.val > no.prob
    else
        Data.Lotto.Type = 2; % group 2 = no.prob > no.val
    end
    val = [5 10 20 40 80];
    prob = [.2 .35 .5 .65 .8];
    comb = combvec(val,prob);
    for n=1:length(comb(1,:))
        repIdx(n,:) = find(Data.Lotto.lottoVal == comb(1,n) & Data.Lotto.lottoProb == comb(2,n)*100);
    end
    % returns a vector of trial index for overlapping choice set:
    Data.Lotto.repIdx = reshape(repIdx,[],1); 
    % feed overlapping choice set into model
    choseLottery = Data.Lotto.choseLotto(Data.Lotto.repIdx);
    refVal = ones(size(choseLottery))*5;
    lottoVal = Data.Lotto.lottoVal(Data.Lotto.repIdx);
    refProb = ones(size(choseLottery));
    lottoProb = Data.Lotto.lottoProb(Data.Lotto.repIdx)/100;
    b0 = [1,1,1]; % initial estimate of alpha and gamma
    fittingResults = fitData_SU(choseLottery,refVal,lottoVal,refProb,lottoProb,b0);
    alpha=fittingResults.b(1);
    gamma=fittingResults.b(2);
    pweight=fittingResults.b(3);
    Data.Lotto.fitRep=[alpha,gamma,pweight];
    % proportion of lottery choices
    Data.Lotto.lottoPropRep = sum(Data.Lotto.choseLotto(Data.Lotto.repIdx))...
        /length(Data.Lotto.choseLotto(Data.Lotto.repIdx));
    %% Feature Bias Index (FBI) 
    % FBI =  (proportion of the trial spent attending to value - proportion
    % of the trial spent attending to probability)
    FBItrial = (Data.Lotto.valDisplayFrames - Data.Lotto.probDisplayFrames)./Data.Lotto.totalFrames;
    Data.Lotto.FBImean = mean(FBItrial);
    %% Save
    save(Data.filename,'Data')
end