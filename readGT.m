function groupList = readGT(vars)
    if length(vars) == 2
        if startsWith(vars{1},'eda')
            oad('output/conditionSchool_EDA_IBI.mat','groupList');
        else
            load('output/conditionSchool_IBI_EDA.mat','groupList');
        end
    else
        if startsWith(vars{1},'eda')
            load('output/conditionSchool_EDA.mat','groupList');
        else
            load('output/conditionSchool_IBI.mat','groupList');
        end
    end
end