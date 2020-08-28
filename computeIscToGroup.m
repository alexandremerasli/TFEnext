function [isc_to_group,accuracy] = computeIscToGroup(isc,groupList,nSubj,nStudiedSubjects,infDiagModes,vars,studiedPeriods,periodsToLookAt,classesToLookAt,varargin)

%    p = inputParser;
%    addOptional(p,'plotFigure',false);
%    parse(p,[isc,groupList,nSubj,nStudiedSubjects,infDiagModes,vars,studiedPeriods,periodsToLookAt,classesToLookAt],varargin{:});

    if nargin > 9
        [varargin{:}] = convertStringsToChars(varargin{:});
    end

    if nargin < 9
        error(message('stats:mdscale:TooFewInputs'));
    end

    paramNames = {'plotFigure'};
    paramDflts = {[]};
    [plotFigure] = internal.stats.parseArgs(paramNames, paramDflts, varargin{:});

    % pre-assignement / parameters
    isc_to_group = cell(1, length(vars)*length(infDiagModes));
    accuracy = zeros(1, length(vars)*length(infDiagModes));

    tmpnSubj = cell(1, length(vars));
    nPeriod = length(periodsToLookAt);
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
            if plotFigure
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
            end
            
            % compute accuracy
            diffWithoutNans = diff(isc_to_group{nComputation}(1 : sum(tmpnSubj{var}(classesToLookAt,periodsToLookAt),'all'),:),1,2); % remove zeros (not studied groups)
            diffWithoutNans = diffWithoutNans(~isnan(diffWithoutNans)); % remove NaNs
            accuracy(nComputation) = sum(diffWithoutNans<0) / length(diffWithoutNans);
        end
    end
end