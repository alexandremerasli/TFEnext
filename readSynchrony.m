function isc = readSynchrony(infDiagModes,vars,studiedPeriods,studiedClasses)
    if length(infDiagModes) == 1
        if infDiagModes
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    load(strcat('schoolInf_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    load(strcat('schoolInf_IBI_EDA',num2str(studiedPeriods),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    load(strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    load(strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        else
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    load(strcat('school_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    load(strcat('school_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    load(strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    load(strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        end
    else
        if infDiagModes(1)
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    load(strcat('schoolInfAvg_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    load(strcat('schoolInfAvg_IBI_EDA',num2str(studiedPeriods),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    load(strcat('schoolInfAvg_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    load(strcat('schoolInfAvg_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        else
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    load(strcat('schoolAvgInf_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    load(strcat('schoolAvgInf_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    load(strcat('schoolAvgInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    load(strcat('schoolAvgInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        end
    end
end