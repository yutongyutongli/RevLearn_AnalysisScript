function [medDiff,UpVal,DownVal,pVal] = permutation(permNum,ListOne,ListTwo,RmedDiff)
    for i = [1:permNum]
    % create permutation list
    DataList = [];
    GroupList = [];
    DataList(:,1) = ListOne;
    GroupList(1:length(ListOne),1)= 1;
    DataList(length(ListOne)+1:100,1) = ListTwo;
    GroupList(length(ListOne)+1:100,1) = 2;
    % randomize group labels
    permGrp = GroupList(randperm(length(GroupList(:,1))));
    % recombine the Data and Group list
    permList(:,1) = DataList; 
    permGrpOne = DataList(permGrp(:,1)==1);
    permGrpTwo = DataList(permGrp(:,1)==2);
    medOne = median(permGrpOne);
    medTwo = median(permGrpTwo);
    medDiff(i) = medOne - medTwo;
    end
    sortPerm = sort(medDiff);
    UpVal = sortPerm(0.025*permNum);
    DownVal = sortPerm(0.975*permNum);
    pVal = sum(abs(RmedDiff)<=abs(medDiff))/permNum;
end