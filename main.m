readISC = true; % keep it false to not overwrite current isc data

% parameters
infDiag = false; % normal synchrony computation or averaging (to avoid inf)
studiedPeriods = [1,2,3,4,5,6,7];
studiedClasses = [1,4];
% studiedClasses = 1:4;

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

% fill data matrix
for var = 1 : length(vars)
    % load data
    path2file = fullfile(filepath, ['for_synchrony_', vars{var}]);
    load(path2file);
    
    % load studied data
    sig_all_studied = sig_firstclass_all(studiedClasses,studiedPeriods);

    % compute number of classes and subjects
    nSubj = cellfun(@(x) length(x), sig_firstclass_all); % number of subjects for each class and each period
    nClass = length(studiedClasses);
    nStudiedSubjects = max(sum(nSubj(studiedClasses,studiedPeriods),1)); % number of studied subjects according to chosen classes and periods
    
    % pre-assign variables for faster processing
    nTime(var) = T*var_srate(var); % 50 min x 60 sec / min x 32 samp / sec %nb of times for each modality
    nPeriod = length(studiedPeriods);
    data{var} = nan(nStudiedSubjects, nTime(var),nPeriod); % will become physiological signals for each chosen modality
    groupList{var} = zeros(nStudiedSubjects, 1); % will become IDs for each chosen modality
    
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
    isc = cell(1, length(vars));
    for var = 1 : length(vars)
        isc{var} = nan(nStudiedSubjects, nStudiedSubjects,nPeriod,nPeriod);
        for n1 = 1 : nStudiedSubjects
            disp(n1/nStudiedSubjects*100); % percentage of done computation
            for period1 = 1:nPeriod
                disp(period1); % can be removed
                for period2 = period1 : nPeriod
                    if period2 == period1
                        bound = n1;
                    else
                        bound = 1;
                    end
                    for n2 = bound : nStudiedSubjects
                        % present data in a struct input needed for 'ps_mwa'
                        dat_subj_1.data = data{var}(n1,t0*var_srate(var)+1:t1*var_srate(var),period1);
                        dat_subj_1.samplerate = var_srate(var);

                        dat_subj_2.data = data{var}(n2,t0*var_srate(var)+1:t1*var_srate(var),period2);
                        dat_subj_2.samplerate = var_srate(var);

                        r = ps_mwa(dat_subj_1, dat_subj_2, ...
                                'CorWindow', 15, ...
                                'CorStep', 1); % synchrony signal

                        if infDiag
                            isc{var}(n1,n2,period1,period2) = r.overall; % usual way to compute synchrony
                        else
                            isc{var}(n1,n2,period1,period2) = nanmean(r.data(8:end-8)); % To avoid Inf on diagonal
                        end
                    end
                end
            end 
        end
    end
end
%% WRITE SYNCHRONY MATRIX IN CSV
if ~readISC % ISC has not been written yet, so we can do it according to 
            % chosen modalities and the way of computing synchrony. 
    if infDiag
        if length(vars) == 2
            if startsWith(vars{1},'eda')
                writematrix(isc{1},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                writematrix(isc{2},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                save(strcat('schoolInf_EDA_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.mat'),'isc');
            else
                writematrix(isc{1},strcat('schoolInf_IBI',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                writematrix(isc{2},strcat('schoolInf_EDA',num2str(studiedPeriods),'_',num2str(studiedClasses),'.csv'));
                save(strcat('schoolInf_IBI_EDA',num2str(studiedPeriods),'.mat'),'isc');
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
end

% ISC matrix has just been stored in csv, we can safely change isc matrix
% or use the previously computed one.
readISC = true;

%%
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
    if infDiag
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
end
%% COMPUTE SYNCHRONY WITH OWN GROUP AND OTHER GROUP

% pre-assignement / parameters
isc_to_group = cell(1, length(vars));
accuracy = zeros(1, length(vars));

% periodsToLookAt = 1; % it is possible to select a subset of studied periods
periodsToLookAt = studiedPeriods(4);
% classesToLookAt = [1,2,4,5,17];
classesToLookAt = studiedClasses;
% classesToLookAt = [1,4]; % it is possible to select a subset of studied classes

tmpnSubj = zeros(size(nSubj)); % number of subjects for each class and each period (can be different from nSubj if classesToLookAt != studiedClasses or periodsToLookAt != studiedPeriods) 
tmpnSubj(classesToLookAt,periodsToLookAt) = nSubj(classesToLookAt,periodsToLookAt);
nPeriod = length(periodsToLookAt);
nClass = length(classesToLookAt);
symIsc = cell(1, length(vars)); % make isc symmetrical to well compute isc_to_group
groupListStudied = cell(1, length(vars)); % ground truth for selected classes

for var = 1 : length(vars)
    groupListStudied{var} = groupList{var}(find(sum(groupList{var}==classesToLookAt,2)));
    
%     for i=1:length(classesToLookAt)
%         classesToLookAt(i) = find(studiedClasses==classesToLookAt(i));
%     end
%     for i=1:length(periodsToLookAt)
%         periodsToLookAt(i) = find(studiedPeriods==periodsToLookAt(i));
%     end
    
    % make isc symmetrical (use a trick to remove NaNs) to well compute isc_to_group
    symIsc{var} = isc{var};
    symIsc{var}(isnan(symIsc{var})) = 1j;
    symIsc{var} = (symIsc{var}+permute(symIsc{var},[2 1 3 4]))/2;    
    NaNsLocation = imag(symIsc{var});
    symIsc{var} = real(symIsc{var});
    symIsc{var}(NaNsLocation==0) = nan;
    symIsc{var} = isc{var};
    
    subj = 0;
    isc_to_group{var} = zeros(nStudiedSubjects*nPeriod, 2);
    for i = 1:nPeriod
        period = periodsToLookAt(i); % studied period
        periodIdx = find(studiedPeriods==periodsToLookAt); % period index in isc
        for n1 = 1 : sum(tmpnSubj(classesToLookAt,period))
            subj = subj + 1;
            % own group
            if find(groupListStudied{var} == groupListStudied{var}(n1)) == n1
                isc_to_group{var}(subj,1) = 1;
            else
                isc_to_group{var}(subj,1) = nanmean(symIsc{var}(n1, setdiff(find(groupListStudied{var} == groupListStudied{var}(n1)), n1), periodIdx, periodIdx)); % ISC within group
            end
            % other group
            if (length(periodsToLookAt) > 1)
                isc_to_group{var}(subj,2) = nanmean(symIsc{var}(n1, groupListStudied{var} ~= groupListStudied{var}(n1),setdiff(1:nPeriod, periodIdx),setdiff(1:nPeriod, periodIdx)),'all'); % ISC between groups
            else
                isc_to_group{var}(subj,2) = nanmean(symIsc{var}(n1, groupListStudied{var} ~= groupListStudied{var}(n1)),'all'); % ISC between groups
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
    for n1 = 1 : length(isc_to_group{var})
        plot([1, 2], isc_to_group{var}(n1, :), ...
                            'Color', color{(diff(isc_to_group{var}(n1,:))>0) +1},... % well classified -> red
                            'Marker', '.', ...
                            'LineWidth', 1.5, ...
                            'LineStyle', linestyle{(diff(isc_to_group{var}(n1,:))>0) +1}, ... % misclassified -> blue
                            'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', color{(diff(isc_to_group{var}(n1,:))>0) +1}, ...
                            'MarkerSize', 15);
    end
    
    % compute accuracy
    diffWithoutNans = diff(isc_to_group{var}(1 : sum(tmpnSubj(classesToLookAt,periodsToLookAt),'all'),:),1,2); % remove zeros (not studied groups)
    diffWithoutNans = diffWithoutNans(~isnan(diffWithoutNans)); % remove NaNs
    accuracy(var) = sum(diffWithoutNans<0) / length(diffWithoutNans);
end

%% COMPUTE SYNCHRONY WITH OWN GROUP AND OTHER GROUP FOR EACH 2-CLASSES COMBINATION

if length(studiedClasses) == 17
    nClassesCombination = 2;
    allClasses = nchoosek(1:17,nClassesCombination); % all nClassesCombination permutations
    for i=1:length(allClasses)
        % pre-assignement / parameters
        isc_to_group = cell(1, length(vars));
        accuracy = zeros(1, length(vars));

        periodsToLookAt = studiedPeriods(1); % Look at only one period
        classesToLookAt = allClasses(i,:);
        tmpnSubj = zeros(size(nSubj));
        tmpnSubj(classesToLookAt,periodsToLookAt) = nSubj(classesToLookAt,periodsToLookAt);
        nPeriod = length(periodsToLookAt);
        nClass = length(classesToLookAt);
        symIsc = cell(1, length(vars));
        groupListStudied = cell(1, length(vars));
        for var = 1 : length(vars)
            groupListStudied{var} = groupList{var}(find(sum(groupList{var}==classesToLookAt,2)));
    
        %     for i=1:length(classesToLookAt)
        %         classesToLookAt(i) = find(studiedClasses==classesToLookAt(i));
        %     end
        %     for i=1:length(periodsToLookAt)
        %         periodsToLookAt(i) = find(studiedPeriods==periodsToLookAt(i));
        %     end

            % make isc symmetrical (use a trick to remove NaNs) to well compute isc_to_group
            symIsc{var} = isc{var};
            symIsc{var}(isnan(symIsc{var})) = 1j;
            symIsc{var} = (symIsc{var}+permute(symIsc{var},[2 1 3 4]))/2;    
            NaNsLocation = imag(symIsc{var});
            symIsc{var} = real(symIsc{var});
            symIsc{var}(NaNsLocation==0) = nan;
            symIsc{var} = isc{var};

            subj = 0;
            isc_to_group{var} = zeros(nStudiedSubjects*nPeriod, 2);
            for periodIdx = 1:nPeriod
                period = periodsToLookAt(periodIdx);
                for n1 = 1 : sum(tmpnSubj(classesToLookAt,period))
                    subj = subj + 1;
                    % own group
                    if find(groupListStudied{var} == groupListStudied{var}(n1)) == n1
                        isc_to_group{var}(subj,1) = 1;
                    else
                        isc_to_group{var}(subj,1) = nanmean(symIsc{var}(n1, setdiff(find(groupListStudied{var} == groupListStudied{var}(n1)), n1), periodIdx, periodIdx)); % ISC within group
                    end
                    % other group
                    if (length(periodsToLookAt) > 1)
                        isc_to_group{var}(subj,2) = nanmean(symIsc{var}(n1, groupListStudied{var} ~= groupListStudied{var}(n1),setdiff(1:nPeriod, periodIdx),setdiff(1:nPeriod, periodIdx)),'all'); % ISC between groups
                    else
                        isc_to_group{var}(subj,2) = nanmean(symIsc{var}(n1, groupListStudied{var} ~= groupListStudied{var}(n1)),'all'); % ISC between groups
                    end
                end
            end

            % compute accuracy
            diffWithoutNans = diff(isc_to_group{var}(1 : sum(tmpnSubj(classesToLookAt,periodsToLookAt),'all'),:),1,2); % remove zeros (not studied groups)
            diffWithoutNans = diffWithoutNans(~isnan(diffWithoutNans)); % remove NaNs
            accuracy(var) = sum(diffWithoutNans<0) / length(diffWithoutNans);
            acc(i,var) = accuracy(var);
        end
    end
end

%%
if length(studiedClasses) == 17
    figure;
    plot(acc);
end