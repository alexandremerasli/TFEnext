addpath('signalFunctions/')
addpath('usefulFunctions/')

readISC = true; % keep it false to not overwrite current isc data

% parameters
infDiagModes = [false]; % normal synchrony computation or averaging (to avoid inf)
% studiedPeriods = [1,2,3];
studiedPeriods = 1:7;
studiedClasses = [1];
% studiedClasses = 1:17;
% vars = {'eda_all_classes', 'hr_all_classes'};
% vars = {'hr_all_classes', 'eda_all_classes'};
% vars = {'hr_all_classes'};
vars = {'eda_all_classes'};
% downsample_eda = 1;
downsample_eda = 8;
downsample_hr = 1;

[data,groupList,nSubj,nTime,nStudiedSubjects,nClass,nPeriod,infDiagModes,vars,t0,t1,var_srate] = initialize(infDiagModes,vars,studiedPeriods,studiedClasses,downsample_eda,downsample_hr);

%% COMPUTE PHYSIOLOGICAL SYNCHRONY

if readISC % we can safely compute new ISC as it has been stored in CSV
             % or user authorized to overwrite it
    readISC = false;
    isc = computeSynchrony(data,nStudiedSubjects,nPeriod,infDiagModes,vars,t0,t1,var_srate);
end


%%
iscPrevious1 = isc;

%% WRITE SYNCHRONY MATRIX IN CSV
if ~readISC % ISC has not been written yet, so we can do it according to 
            % chosen modalities and the way of computing synchrony. 
    writeSynchrony(isc,infDiagModes,vars,studiedPeriods,studiedClasses);
end

% ISC matrix has just been stored in csv, we can safely change isc matrix
% or use the previously computed one.
readISC = true;

%% Just in case
iscPrevious = isc;

%% WRITE TRUE GROUPS in CSV
if nClass == 17 % write only if every classes are selected to match with Python then
    writeGT(groupList,vars);
end

%% READ SYNCHRONY MATRIX IN CSV 
if readISC % we can overwrite isc matrix if user authorized it, or overwrite 
           % by previously computed matrix according to chosen modalities 
           % and the way of computing synchrony
    isc = readSynchrony(infDiagModes,vars,studiedPeriods,studiedClasses);
    groupList = readGT(vars);
    %groupList = groupList(find(sum(groupList==studiedClasses,2)));
end
%% COMPUTE SYNCHRONY WITH OWN GROUP AND OTHER GROUP
% parameters
% periodsToLookAt = 1; % it is possible to select a subset of studied periods
periodsToLookAt = studiedPeriods;
% classesToLookAt = [1,4]; % it is possible to select a subset of studied classes
classesToLookAt = studiedClasses;

% computation
[isc_to_group,accuracy] = computeIscToGroup(isc,groupList,nSubj,nStudiedSubjects,infDiagModes,vars,studiedPeriods,periodsToLookAt,classesToLookAt,'plotFigure',true);
%% COMPUTE SYNCHRONY WITH OWN GROUP AND OTHER GROUP FOR EACH 2-CLASSES COMBINATION

if length(studiedClasses) == 17
    nClassesCombination = 2;
    allClasses = nchoosek(1:17,nClassesCombination); % all nClassesCombination-classes permutations
    accClasses = zeros(length(allClasses), length(vars)*length(infDiagModes));
    nComputation = 0; % counter just for indexing

    for var = 1 : length(vars)
        for infDiag = infDiagModes
            nComputation = nComputation + 1;
            for i=1:length(allClasses)
                classesToLookAt = allClasses(i,:);
                [isc_to_group,accClasses(i,nComputation)] = computeIscToGroup(isc,groupList,nSubj,nStudiedSubjects,infDiagModes,vars,studiedPeriods,periodsToLookAt,classesToLookAt);
            end
        end
    end
end

%% PLOT ACCURACY FOR EACH COMBINATION OF CLASSES
if length(studiedClasses) == 17
    figure;
    plot(accClasses);
end

%% STUDY PERFORMANCE BETWEEN PERIODS IN ONE GROUP

infDiagModes = false;
vars = {'eda_all_classes'};
studiedPeriods = 1:7;
studiedClasses = studiedClasses(1);

if readISC % we can overwrite isc matrix if user authorized it, or overwrite 
           % by previously computed matrix according to chosen modalities 
           % and the way of computing synchrony
    readSynchrony(infDiagModes,vars,studiedPeriods,studiedClasses);
end

classesToLookAt = studiedClasses;
nPeriodsCombination = 2;
allPeriods = nchoosek(1:6,nPeriodsCombination); % all nPeriodsCombination-periods permutations
accPeriods = zeros(length(allPeriods), length(vars)*length(infDiagModes));
nComputation = 0; % counter just for indexing

for var = 1 : length(vars)
    for infDiag = infDiagModes
        nComputation = nComputation + 1;
        for i=1:length(allPeriods)
            periodsToLookAt = allPeriods(i,:);
            [isc_to_group,accPeriods(i,nComputation)] = computeIscToGroup(isc,groupList,nSubj,nStudiedSubjects,infDiagModes,vars,studiedPeriods,periodsToLookAt,classesToLookAt);
        end
    end
end

%% PLOT ACCURACY FOR EACH COMBINATION OF PERIODS
n = length(studiedPeriods);
if n == 7
    n = 6;
    figure;
    legends = string(zeros(1,n-1));
    for i = 1:n-2
        hold on;
        plot(i+1:n,accPeriods(((i-1)*n - i*(i-1)/2 + 1):(i*n - i*(i+1)/2)),'LineWidth',2);
        legends(i) = strcat("period ",num2str(i));
    end
    i = n - 1;
    scatter(i+1:n,accPeriods(((i-1)*n - i*(i-1)/2 + 1):(i*n - i*(i+1)/2)),'LineWidth',2,'Marker','x');
    legends(i) = strcat("period ",num2str(i));
    legend(legends);
    xlabel('other period');
    ylabel('accuracy');
end