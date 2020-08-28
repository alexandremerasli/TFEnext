function writeGT(groupList,vars)
    if length(vars) == 2
        if startsWith(vars{1},'eda')
            writematrix(groupList{1},'output/conditionSchool_EDA.csv');
            writematrix(groupList{2},'output/conditionSchool_IBI.csv');
            save('output/conditionSchool_EDA_IBI.mat','groupList');
        else
            writematrix(groupList{1},'output/conditionSchool_IBI.csv');
            writematrix(groupList{2},'output/conditionSchool_EDA.csv');
            save('output/conditionSchool_IBI_EDA.mat','groupList');
        end
    else
        if startsWith(vars{1},'eda')
            writematrix(groupList{1},'output/conditionSchool_EDA.csv');
            save('output/conditionSchool_EDA.mat','groupList');
        else
            writematrix(groupList{1},'output/conditionSchool_IBI.csv');
            save('output/conditionSchool_IBI.mat','groupList');
        end
    end
end