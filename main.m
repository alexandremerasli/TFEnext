readISC = true; % keep it false to not overwrite current isc data

% parameters
infDiagModes = [false]; % normal synchrony computation or averaging (to avoid inf)
% studiedPeriods = [2];
studiedPeriods = 1:7;
% studiedClasses = [1,4];
studiedClasses = 1:17;

% choose modality (do not run both for now, as different labels, will
% produce errors)
% vars = {'eda_all_classes', 'hr_all_classes'};
% vars = {'hr_all_classes', 'eda_all_classes'};
% vars = {'hr_all_classes'};
vars = {'eda_all_classes'};

% signal information according to chosen modalities
if length(vars) == 2
    if startsWith(vars{1},'eda')
        var_srate = [32,1];
    else
        var_srate = [1,32];
    end
else
    if startsWith(vars{1},'eda')
        var_srate = [32];
    else
        var_srate = [1];
    end
end

T = 50 * 60; % recording time [s]
t0 = 5 * 60; % start time of analysis [s];
t1 = (50 - 5) * 60; % end time of analysis [s]

filepath = './data';

% pre-assign variables for faster processing
nTime = zeros(1, length(vars));
data = cell(1, length(vars)); % physiological signals
groupList = cell(1, length(vars)); % IDs

nSubj = cell(1, length(vars));
nClass = zeros(1, length(vars));
nStudiedSubjects = zeros(1, length(vars));

% fill data matrix
for var = 1 : length(vars)
    % load data
    path2file = fullfile(filepath, ['for_synchrony_', vars{var}]);
    load(path2file);
    
    % load studied data
    sig_all_studied = sig_firstclass_all(studiedClasses,studiedPeriods);

    % compute number of classes and subjects
    nSubj{var} = cellfun(@(x) length(x), sig_firstclass_all); % number of subjects for each class and each period
    nClass(var) = length(studiedClasses);
    nStudiedSubjects(var) = max(sum(nSubj{var}(studiedClasses,studiedPeriods),1)); % number of studied subjects according to chosen classes and periods
    
    % pre-assign variables for faster processing
    nTime(var) = T*var_srate(var); % 50 min x 60 sec / min x 32 samp / sec %nb of times for each modality
    nPeriod = length(studiedPeriods);
    data{var} = nan(nStudiedSubjects(var), nTime(var),nPeriod); % will become physiological signals for each chosen modality
    groupList{var} = zeros(nStudiedSubjects(var), 1); % will become IDs for each chosen modality
    
    % fill matrix 'data' with responses
    subj = 0;
    for cl = 1 : nClass
        for subj_in_cl = 1 : length(sig_firstclass_all{studiedClasses(cl)})
            subj = subj + 1;
            groupList{var}(subj) = studiedClasses(cl);
            for period=1:nPeriod
                if isempty(sig_all_studied{cl,period})
                    continue
                elseif length(sig_all_studied{cl,period}{subj_in_cl}) ~= nTime(var) % skip data with different number of times
                    continue
                else
                    data{var}(subj, :,period) = sig_all_studied{cl,period}{subj_in_cl}; % fill data
                end
            end
        end 
    end
end

%% COMPUTE PHYSIOLOGICAL SYNCHRONY

if readISC % we can safely compute new ISC as it has been stored in CSV
             % or user authorized to overwrite it
    readISC = false;
    nComputation = 0; % counter just for indexing
    isc = cell(1, length(vars)*length(infDiagModes));
    for var = 1 : length(vars)
        disp(vars{var})
        for infDiag = infDiagModes
            disp(infDiag)
            nComputation = nComputation + 1; 
            isc{nComputation} = nan(nStudiedSubjects(var), nStudiedSubjects(var),nPeriod,nPeriod);
            for n1 = 1 : nStudiedSubjects(var)
                disp(n1/nStudiedSubjects(var)*100); % percentage of done computation
                for period1 = 1:nPeriod
                    disp(period1); % can be removed
                    for period2 = period1 : nPeriod
                        if period2 == period1
                            bound = n1;
                        else
                            bound = 1;
                        end
                        for n2 = bound : nStudiedSubjects(var)
                            % present data in a struct input needed for 'ps_mwa'
                            dat_subj_1.data = data{var}(n1,t0*var_srate(var)+1:t1*var_srate(var),period1);
                            dat_subj_1.samplerate = var_srate(var);

                            dat_subj_2.data = data{var}(n2,t0*var_srate(var)+1:t1*var_srate(var),period2);
                            dat_subj_2.samplerate = var_srate(var);

                            r = ps_mwa(dat_subj_1, dat_subj_2, ...
                                    'CorWindow', 15, ...
                                    'CorStep', 1); % synchrony signal

                            if infDiag
                                isc{nComputation}(n1,n2,period1,period2) = r.overall; % usual way to compute synchrony
                            else
                                isc{nComputation}(n1,n2,period1,period2) = nanmean(r.data(8:end-8)); % To avoid Inf on diagonal
                            end
                        end
                    end
                end 
            end
        end
    end
end


%%
iscPrevious1 = isc;

%% WRITE SYNCHRONY MATRIX IN CSV
if ~readISC % ISC has not been written yet, so we can do it according to 
            % chosen modalities and the way of computing synchrony. 
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
                    save(strcat('school_EDA_IBI',num2str(studiedPeriods),'.mat'),'_',num2str(studiedClasses),'isc');
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

% ISC matrix has just been stored in csv, we can safely change isc matrix
% or use the previously computed one.
readISC = true;

%% Just in case
iscPrevious = isc;

%% WRITE TRUE GROUPS in CSV
if nClass == 17 % write only if every classes are selected to match with Python then
    if length(vars) == 2
        if startsWith(vars{1},'eda')
            writematrix(groupList{1},'conditionSchool_EDA.csv');
            writematrix(groupList{2},'conditionSchool_IBI.csv');
        else
            writematrix(groupList{1},'conditionSchool_IBI.csv');
            writematrix(groupList{2},'conditionSchool_EDA.csv');
        end
    else
        if startsWith(vars{1},'eda')
            writematrix(groupList{1},'conditionSchool_EDA.csv');
        else
            writematrix(groupList{1},'conditionSchool_IBI.csv');
        end
    end
end

%% READ SYNCHRONY MATRIX IN CSV 
if readISC % we can overwrite isc matrix if user authorized it, or overwrite 
           % by previously computed matrix according to chosen modalities 
           % and the way of computing synchrony
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
%% COMPUTE SYNCHRONY WITH OWN GROUP AND OTHER GROUP

% pre-assignement / parameters
isc_to_group = cell(1, length(vars)*length(infDiagModes));
accuracy = zeros(1, length(vars)*length(infDiagModes));

% periodsToLookAt = 1; % it is possible to select a subset of studied periods
periodsToLookAt = studiedPeriods;
% classesToLookAt = [1,4]; % it is possible to select a subset of studied classes
classesToLookAt = studiedClasses;

tmpnSubj = cell(1, length(vars));
nPeriod = length(periodsToLookAt);
nClass = length(classesToLookAt);
symIsc = cell(1, length(vars)*length(infDiagModes)); % make isc symmetrical to well compute isc_to_group
groupListStudied = cell(1, length(vars)); % ground truth for selected classes

nComputation = 0; % counter just for indexing

for var = 1 : length(vars)
    for infDiag = infDiagModes
        nComputation = nComputation + 1;
        tmpnSubj{var} = zeros(size(nSubj)); % number of subjects for each class and each period (can be different from nSubj if classesToLookAt != studiedClasses or periodsToLookAt != studiedPeriods) 
        tmpnSubj{var}(classesToLookAt,periodsToLookAt) = nSubj{var}(classesToLookAt,periodsToLookAt);

        groupListStudied{var} = groupList{var}(find(sum(groupList{var}==classesToLookAt,2)));

        % make isc symmetrical (use a trick with NaNs) to well compute isc_to_group
        symIsc{nComputation} = isc{nComputation};
        symIsc{nComputation}(isnan(symIsc{nComputation})) = 1j;
        symIsc{nComputation} = (symIsc{nComputation}+permute(symIsc{nComputation},[2 1 3 4]));    
        NaNsLocation = imag(symIsc{nComputation});
        symIsc{nComputation} = real(symIsc{nComputation});
        symIsc{nComputation}(NaNsLocation==2) = nan;

        subj = 0;
        isc_to_group{nComputation} = zeros(nStudiedSubjects(var)*nPeriod, 2);
        for p = 1:nPeriod
            period = periodsToLookAt(p); % studied period
            periodIdx = find(studiedPeriods==period); % period index in isc
            for n1 = 1 : sum(tmpnSubj{var}(classesToLookAt,period))
                subj = subj + 1;
                % own group
                if find(groupListStudied{var} == groupListStudied{var}(n1)) == n1
                    isc_to_group{nComputation}(subj,1) = 1;
                else
                    isc_to_group{nComputation}(subj,1) = nanmean(symIsc{nComputation}(n1, setdiff(find(groupListStudied{var} == groupListStudied{var}(n1)), n1), periodIdx, periodIdx)); % ISC within group
                end
                % other group
                if (length(periodsToLookAt) > 1)
                    isc_to_group{nComputation}(subj,2) = nanmean(symIsc{nComputation}(n1, groupListStudied{var} ~= groupListStudied{var}(n1),setdiff(1:nPeriod, periodIdx),setdiff(1:nPeriod, periodIdx)),'all'); % ISC between groups
                else
                    isc_to_group{nComputation}(subj,2) = nanmean(symIsc{nComputation}(n1, groupListStudied{var} ~= groupListStudied{var}(n1),periodIdx,periodIdx),'all'); % ISC between groups
                end
            end
        end

        % figure
        color = {'Red', 'Blue'};
        linestyle = {'-', '--'};

        fig = figure;
        ax = axes('Parent', fig, ...
            'xlim', [0, 3]);
        hold on;
        for n1 = 1 : length(isc_to_group{nComputation})
            plot([1, 2], isc_to_group{nComputation}(n1, :), ...
                                'Color', color{(diff(isc_to_group{nComputation}(n1,:))>0) +1},... % well classified -> red
                                'Marker', '.', ...
                                'LineWidth', 1.5, ...
                                'LineStyle', linestyle{(diff(isc_to_group{nComputation}(n1,:))>0) +1}, ... % misclassified -> blue
                                'MarkerFaceColor', 'none', ...
                                'MarkerEdgeColor', color{(diff(isc_to_group{nComputation}(n1,:))>0) +1}, ...
                                'MarkerSize', 15);
        end

        % compute accuracy
        diffWithoutNans = diff(isc_to_group{nComputation}(1 : sum(tmpnSubj{var}(classesToLookAt,periodsToLookAt),'all'),:),1,2); % remove zeros (not studied groups)
        diffWithoutNans = diffWithoutNans(~isnan(diffWithoutNans)); % remove NaNs
        accuracy(nComputation) = sum(diffWithoutNans<0) / length(diffWithoutNans);
    end
end
%% COMPUTE SYNCHRONY WITH OWN GROUP AND OTHER GROUP FOR EACH 2-CLASSES COMBINATION

if length(studiedClasses) == 17
    nClassesCombination = 2;
    allClasses = nchoosek(1:17,nClassesCombination); % all nClassesCombination-classes permutations
    acc = zeros(length(allClasses), length(vars)*length(infDiagModes));
    for i=1:length(allClasses)
        % pre-assignement / parameters
        isc_to_group = cell(1, length(vars)*length(infDiagModes));
        accuracy = zeros(1, length(vars)*length(infDiagModes));
        periodsToLookAt = studiedPeriods(1); % Look at only one period
        classesToLookAt = allClasses(i,:);
        tmpnSubj = cell(1, length(vars));
        nPeriod = length(periodsToLookAt);
        nClass = length(classesToLookAt);
        symIsc = cell(1, length(vars)*length(infDiagModes));
        groupListStudied = cell(1, length(vars));
        nComputation = 0; % counter just for indexing

        for var = 1 : length(vars)
            for infDiag = infDiagModes
                nComputation = nComputation + 1;
                tmpnSubj{var} = zeros(size(nSubj)); % number of subjects for each class and each period (can be different from nSubj if classesToLookAt != studiedClasses or periodsToLookAt != studiedPeriods) 
                tmpnSubj{var}(classesToLookAt,periodsToLookAt) = nSubj{var}(classesToLookAt,periodsToLookAt);

                groupListStudied{var} = groupList{var}(find(sum(groupList{var}==classesToLookAt,2)));

                % make isc symmetrical (use a trick with NaNs) to well compute isc_to_group
                symIsc{nComputation} = isc{nComputation};
                symIsc{nComputation}(isnan(symIsc{nComputation})) = 1j;
                symIsc{nComputation} = (symIsc{nComputation}+permute(symIsc{nComputation},[2 1 3 4]));    
                NaNsLocation = imag(symIsc{nComputation});
                symIsc{nComputation} = real(symIsc{nComputation});
                symIsc{nComputation}(NaNsLocation==2) = nan;

                subj = 0;
                isc_to_group{nComputation} = zeros(nStudiedSubjects(var)*nPeriod, 2);
                for p = 1:nPeriod
                    period = periodsToLookAt(p); % studied period
                    periodIdx = find(studiedPeriods==period); % period index in isc
                    for n1 = 1 : sum(tmpnSubj{var}(classesToLookAt,period))
                        subj = subj + 1;
                        % own group
                        if find(groupListStudied{var} == groupListStudied{var}(n1)) == n1
                            isc_to_group{nComputation}(subj,1) = 1;
                        else
                            isc_to_group{nComputation}(subj,1) = nanmean(symIsc{nComputation}(n1, setdiff(find(groupListStudied{var} == groupListStudied{var}(n1)), n1), periodIdx, periodIdx)); % ISC within group
                        end
                        % other group
                        if (length(periodsToLookAt) > 1)
                            isc_to_group{nComputation}(subj,2) = nanmean(symIsc{nComputation}(n1, groupListStudied{var} ~= groupListStudied{var}(n1),setdiff(1:nPeriod, periodIdx),setdiff(1:nPeriod, periodIdx)),'all'); % ISC between groups
                        else
                            isc_to_group{nComputation}(subj,2) = nanmean(symIsc{nComputation}(n1, groupListStudied{var} ~= groupListStudied{var}(n1),periodIdx,periodIdx),'all'); % ISC between groups
                        end
                    end
                end
                
                % compute accuracy
                diffWithoutNans = diff(isc_to_group{nComputation}(1 : sum(tmpnSubj{var}(classesToLookAt,periodsToLookAt),'all'),:),1,2); % remove zeros (not studied groups)
                diffWithoutNans = diffWithoutNans(~isnan(diffWithoutNans)); % remove NaNs
                accuracy(nComputation) = sum(diffWithoutNans<0) / length(diffWithoutNans);
                acc(i,nComputation) = accuracy(nComputation);
            end
        end
    end
end

%% PLOT ACCURACY FOR EACH COMBINATION OF CLASSES
if length(studiedClasses) == 17
    figure;
    plot(acc);
end