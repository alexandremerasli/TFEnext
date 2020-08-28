function isc = computeSynchrony(data,nStudiedSubjects,nPeriod,infDiagModes,vars,t0,t1,var_srate)
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