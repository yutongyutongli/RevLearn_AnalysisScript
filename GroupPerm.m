function [UpVal,DownVal] = GroupPerm(permNum,ListOne)
    for i = [1:permNum]
    % create permutation list
    DataList = [];
    DataList(:,1) = ListOne;
    y = datasample(ListOne,size(ListOne,2));
    medOne(i) = median(y);
    end
    sortPerm = sort(medOne);
    UpVal = sortPerm(0.025*permNum);
    DownVal = sortPerm(0.975*permNum);
end