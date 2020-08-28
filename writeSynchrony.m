function writeSynchrony(isc,infDiagModes,vars,studiedPeriods,studiedClasses)
    if length(infDiagModes) == 1
        if infDiagModes
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolInf_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolInf_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        else
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('school_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('school_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        end
    else
        if infDiagModes(1)
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{3},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{4},strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolInfAvg_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{3},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{4},strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolInfAvg_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolInfAvg_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolInfAvg_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        else
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{3},strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{4},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolAvgInf_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{3},strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{4},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolAvgInf_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolAvgInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('schoolAvgInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        end
    end
    end