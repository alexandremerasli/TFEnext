function groupList = readGT(vars)
    if length(vars) == 2
        if startsWith(vars{1},'eda')
            oad('conditionSchool_EDA_IBI.mat','groupList');
        else
            load('conditionSchool_IBI_EDA.mat','groupList');
        end
    else
        if startsWith(vars{1},'eda')
            load('conditionSchool_EDA.mat','groupList');
        else
            load('conditionSchool_IBI.mat','groupList');
        end
    end
end