function [data,groupList,nSubj,nTime,nStudiedSubjects,nClass,nPeriod,infDiagModes,vars,t0,t1,var_srate] = initialize(infDiagModes,vars,studiedPeriods,studiedClasses,downsample_eda,downsample_hr)
    % signal information according to chosen modalities
    eda_freq = 32 / downsample_eda;
    hr_freq = 1 / downsample_hr;

    if length(vars) == 2
        if startsWith(vars{1},'eda')
            var_srate = [eda_freq,hr_freq];
        else
            var_srate = [hr_freq,eda_freq];
        end
    else
        if startsWith(vars{1},'eda')
            var_srate = [eda_freq];
        else
            var_srate = [hr_freq];
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

        % downsampling factor
        if startsWith(vars(var),'eda')
            downsample_var = downsample_eda;
        else
            downsample_var = downsample_hr;
        end

        % fill matrix 'data' with responses
        subj = 0;
        for cl = 1 : nClass
            for subj_in_cl = 1 : length(sig_firstclass_all{studiedClasses(cl)})
                subj = subj + 1;
                groupList{var}(subj) = studiedClasses(cl);
                for period=1:nPeriod
                    if isempty(sig_all_studied{cl,period})
                        continue
                    else
                        if isempty(sig_all_studied{cl,period}{subj_in_cl})
                            continue
                        else
                            sig = downsample(sig_all_studied{cl,period}{subj_in_cl},downsample_var);
                            if length(sig) ~= nTime(var) % skip data with different number of times
                                continue
                            else
                                data{var}(subj, :,period) = sig; % fill data
                            end
                        end
                    end
                end
            end 
        end
    end
end