readISC = true; % Keep it false to not overwrite current isc data
infDiag = false;
periodStudied = 1:2;

% vars = {'eda_all_classes', 'hr_all_classes'};
% vars = {'hr_all_classes', 'eda_all_classes'};
% vars = {'hr_all_classes'};
vars = {'eda_all_classes'};

% var_srate = [32, 1];
% var_srate = [1,32]; 
% var_srate = [1];
var_srate = [32];

T = 50 * 60; % recording time [s]
t0 = 5 * 60; % start time of analysis [s];
t1 = (50 - 5) * 60; % end time of analysis [s]

filepath = './data';

% pre-assign variables for faster processing
nPeriod = zeros(1, length(vars));
for var = 1 : length(vars)
    path2file = fullfile(filepath, ['for_synchrony_', vars{var}]);
    load(path2file);
    sig_firstclass_all = sig_firstclass_all(:,periodStudied);
    nPeriod = size(sig_firstclass_all,2);
end
nSubj = zeros(nPeriod, length(vars));
nTime = zeros(1, length(vars));
data = cell(1, length(vars));
groupList = cell(1, length(vars));

for var = 1 : length(vars)
    % load data
    path2file = fullfile(filepath, ['for_synchrony_', vars{var}]);
    load(path2file);
    sig_firstclass_all = sig_firstclass_all(:,periodStudied);
    
    % sort participants
    %[nPeriod,idx] = sort(sum(cellfun(@(x) length(x), sig_firstclass_all)>0,2),'descend');
    
    
    % compute number of subjects
    nSubj(:,var) = sum(cellfun(@(x) length(x), sig_firstclass_all)); %nb of subjects for each modality
    nTime(var) = T*var_srate(var); % 50 min x 60 sec / min x 32 samp / sec %nb of times for each modality
    data{var} = nan(nSubj(1,var), nTime(var),nPeriod); % will become physiological signals

    nSubjectStudied = 12;
    nClusters = cellfun(@(x) length(x), sig_firstclass_all);
    nClusters = cumsum(nClusters(:,1));
    nClusters = find(nClusters==nSubjectStudied);
    % nSubjectStudied = nSubj(1,var);

    groupList{var} = zeros(nSubjectStudied, 1); % will become IDs
    
    % fill matrix 'data' with responses
    subj = 0;
    
    for cl = 1 : nClusters
        for subj_in_cl = 1 : length(sig_firstclass_all{cl})
            subj = subj + 1;
            groupList{var}(subj) = cl;
            for period=1:nPeriod
                if isempty(sig_firstclass_all{cl,period})
                    continue
                elseif length(sig_firstclass_all{cl,period}{subj_in_cl}) ~= nTime(var) % skip data with different number of times
                    %nLessSubj(period,var) = nLessSubj(period,var) - 1;
                    continue
                else
                    data{var}(subj, :,period) = sig_firstclass_all{cl,period}{subj_in_cl}; % fill data
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
        isc{var} = nan(nSubjectStudied, nSubjectStudied,nPeriod,nPeriod);
        for n1 = 1 : nSubjectStudied
            disp(n1/nSubjectStudied*100)
            for period1 = 1 : nPeriod
                period1
                for period2 = period1 : nPeriod
                    if period2 == period1
                        bound = n1;
                    else
                        bound = 1;
                    end
                    for n2 = bound : nSubjectStudied
                    % disp(n2)
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
if ~readISC % ISC has not been written yet, so we can do it. 
    if isc{1}(1,1,1) == Inf
        if length(vars) == 2
            if startsWith(vars{1},'eda')
                writematrix(isc{1},strcat('schoolInf_EDA',num2str(nSubjectStudied),'.csv'));
                writematrix(isc{2},strcat('schoolInf_IBI',num2str(nSubjectStudied),'.csv'));
                save(strcat('schoolInf_EDA_IBI',num2str(nSubjectStudied),'.mat'),'isc');
            else
                writematrix(isc{1},strcat('schoolInf_IBI',num2str(nSubjectStudied),'.csv'));
                writematrix(isc{2},strcat('schoolInf_EDA',num2str(nSubjectStudied),'.csv'));
                save(strcat('schoolInf_IBI_EDA',num2str(nSubjectStudied),'.mat'),'isc');
            end
        else
            if startsWith(vars{1},'eda')
                writematrix(isc{1},strcat('schoolInf_EDA',num2str(nSubjectStudied),'.csv'));
                save(strcat('schoolInf_EDA',num2str(nSubjectStudied),'.mat'),'isc');
            else
                writematrix(isc{1},strcat('schoolInf_IBI',num2str(nSubjectStudied),'.csv'));
                save(strcat('schoolInf_IBI',num2str(nSubjectStudied),'.mat'),'isc');
            end
        end
    else
        if length(vars) == 2
            if startsWith(vars{1},'eda')
                writematrix(isc{1},strcat('school_EDA',num2str(nSubjectStudied),'.csv'));
                writematrix(isc{2},strcat('school_IBI',num2str(nSubjectStudied),'.csv'));
                save(strcat('school_EDA_IBI',num2str(nSubjectStudied),'.mat'),'isc');
            else
                writematrix(isc{1},strcat('school_IBI',num2str(nSubjectStudied),'.csv'));
                writematrix(isc{2},strcat('school_EDA',num2str(nSubjectStudied),'.csv'));
                save(strcat('school_IBI_EDA',num2str(nSubjectStudied),'.mat'),'isc');
            end
        else
            if startsWith(vars{1},'eda')
                writematrix(isc{1},strcat('school_EDA',num2str(nSubjectStudied),'.csv'));
                save(strcat('school_EDA',num2str(nSubjectStudied),'.mat'),'isc');
            else
                writematrix(isc{1},strcat('school_IBI',num2str(nSubjectStudied),'.csv'));
                save(strcat('school_IBI',num2str(nSubjectStudied),'.mat'),'isc');
            end
        end
    end
end
% ISC matrix has just been stored in csv, we can safely change isc matrix
% or use the previously computed one.
readISC = true;

%% WRITE TRUE GROUPS in CSV
% writematrix(groupList{1},'conditionSchool_EDA.csv')
% writematrix(groupList{2},'conditionSchool_IBI.csv')

%% READ SYNCHRONY MATRIX IN CSV 
if readISC
    if infDiag
        if length(vars) == 2
            if startsWith(vars{1},'eda')
                load(strcat('schoolInf_EDA_IBI',num2str(nSubjectStudied),'.mat'),'isc');
            else
                load(strcat('schoolInf_IBI_EDA',num2str(nSubjectStudied),'.mat'),'isc');
            end
        else
            if startsWith(vars{1},'eda')
                load(strcat('schoolInf_EDA',num2str(nSubjectStudied),'.mat'),'isc');
            else
                load(strcat('schoolInf_IBI',num2str(nSubjectStudied),'.mat'),'isc');
            end
        end
    else
        if length(vars) == 2
            if startsWith(vars{1},'eda')
                load(strcat('school_EDA_IBI',num2str(nSubjectStudied),'.mat'),'isc');
            else
                load(strcat('school_IBI_EDA',num2str(nSubjectStudied),'.mat'),'isc');
            end
        else
            if startsWith(vars{1},'eda')
                load(strcat('school_EDA',num2str(nSubjectStudied),'.mat'),'isc');
            else
                load(strcat('school_IBI',num2str(nSubjectStudied),'.mat'),'isc');
            end
        end
    end
end
%% COMPUTE SYNCHRONY WITH OWN GROUP AND OTHER GROUP

isc_to_group = cell(1, length(vars));
accuracy = zeros(1, length(vars));

periodToLookAt = 1;
periodToLookAt = periodStudied;

tmpnSubj = nSubj;
tmpnSubj = zeros(size(tmpnSubj));
tmpnSubj(periodToLookAt,:) = nSubj(periodToLookAt,:);
nPeriod = length(periodToLookAt);

for var = 1 : length(vars)
    %var
    isc_to_group{var} = zeros(sum(tmpnSubj(:,var)), 2);
    for period = periodToLookAt
        %period
        for n1 = 1 : nSubjectStudied
            %n1
            %sum(tmpnSubj(1:period2-1,var))+n1
            isc{var}(n1, setdiff(find(groupList{var} == groupList{var}(n1)), n1), period,period) % disp

            % own group
            isc_to_group{var}(sum(tmpnSubj(1:period-1,var))+n1,1) = nanmean(isc{var}(n1, setdiff(find(groupList{var} == groupList{var}(n1)), n1), period, period)); % ISC within group

            % other group
            if (length(periodToLookAt) > 1)
                isc_to_group{var}(sum(tmpnSubj(1:period-1,var))+n1,2) = nanmean(isc{var}(n1, groupList{var} ~= groupList{var}(n1),setdiff(1:nPeriod, period),setdiff(1:nPeriod, period)),'all'); % ISC between groups
            else
                isc_to_group{var}(sum(tmpnSubj(1:period-1,var))+n1,2) = nanmean(isc{var}(n1, groupList{var} ~= groupList{var}(n1)),'all'); % ISC between groups
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
    for n1 = 1 : sum(tmpnSubj(:,var))
        plot([1, 2], isc_to_group{var}(n1, :), ...
                            'Color', color{(diff(isc_to_group{var}(n1,:))>0) +1},... % well classified -> red
                            'Marker', '.', ...
                            'LineWidth', 1.5, ...
                            'LineStyle', linestyle{(diff(isc_to_group{var}(n1,:))>0) +1}, ... % misclassified -> blue
                            'MarkerFaceColor', 'none', ...
                            'MarkerEdgeColor', color{(diff(isc_to_group{var}(n1,:))>0) +1}, ...
                            'MarkerSize', 15);
    end
    diffWithoutNans = diff(isc_to_group{var}(1 : sum(tmpnSubj(:,var)),:),1,2);
    diffWithoutNans = diffWithoutNans(~isnan(diffWithoutNans));
    accuracy(var) = sum(diffWithoutNans<0) / length(diffWithoutNans);
end