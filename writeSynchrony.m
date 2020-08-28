function writeSynchrony(isc,infDiagModes,vars,studiedPeriods,studiedClasses)
    if length(infDiagModes) == 1
        if infDiagModes
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolInf_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolInf_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        else
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/school_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/school_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        end
    else
        if infDiagModes(1)
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{3},strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{4},strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolInfAvg_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{3},strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{4},strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolInfAvg_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolInfAvg_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolInfAvg_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        else
            if length(vars) == 2
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{3},strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{4},strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolAvgInf_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{3},strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{4},strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolAvgInf_IBI_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            else
                if startsWith(vars{1},'eda')
                    writematrix(isc{1},strcat('output/school_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolAvgInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                else
                    writematrix(isc{1},strcat('output/school_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    writematrix(isc{2},strcat('output/schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                    save(strcat('output/schoolAvgInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
                end
            end
        end
    end
    end