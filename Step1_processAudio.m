%% PROCESS AUDIO
% This script generates speaker level features
%
% Author: Tuka Alhanai CSAIL MIT 2017 Copyright

clear all

% Step 0. You have to have a file that defines the segments for a given speaker
% (this is an invented speaker alignment)
rttmFile = '1234.rttm.txt';


% Step 1. Feature extracted using opensmile - but you can feed in what you like.
show = '1234.lld.csv';
name = '1234';

% init
DATA = table;

for speechCutoff = [0.1]
            
    % LOOP THROUGH ALL SUBJECTS
    for i = 1
        
        i
        
        count = 1;
        
        % import data
        data = readtable(char(show));
        names = fieldnames(data);

        % toss out audspec
        toss = contains(names,'audspec');
        names(toss) = [];

        % last 3 paramters are extra stuff from table definition
        data(:,toss(1:end-3)) = [];
        
        % Step 2. Zscore the following parameters - across the complete file.
        for k = 1:14
            eval(['data.mfcc_sma_',   num2str(k),'_ = zscore(data.mfcc_sma_',    num2str(k),'_);']);
            eval(['data.mfcc_sma_de_',num2str(k),'_ = zscore(data.mfcc_sma_de_', num2str(k),'_);']);
        end
        
        data.pcm_zcr_sma    = zscore(data.pcm_zcr_sma);
        data.pcm_zcr_sma_de = zscore(data.pcm_zcr_sma_de);
        
        data.pcm_RMSenergy_sma    = zscore(data.pcm_RMSenergy_sma);
        data.pcm_RMSenergy_sma_de = zscore(data.pcm_RMSenergy_sma_de);

        data.logHNR_sma    = zscore(data.logHNR_sma);
        data.logHNR_sma_de = zscore(data.logHNR_sma_de);
        
        % Step 3. Smooth the voicing signal.
        dataVoice = envelope(zscore(data.voicingFinalUnclipped_sma),10,'peak');
        
        % import diarized segment file
        rttm = readtable(rttmFile);
        
        % we only care for the segments the subject spoke, as denoted by '-P'.
        ind = contains(cellstr(rttm.Var8),'-P');
        rttm = rttm(ind,:);
        
        % INIT ARRAY
        % this will contain all the speech frames for the subject.
        subsetData_spkr0 = [];

        % LOOP THROUGH EACH SUBJECT INDEX
        % -------------------------------
        for j = 1:height(rttm)

            % GET THE START/END TIMES OF THIS SEGMENT
            startTime = rttm.Var4(j);
            endTime   = rttm.Var4(j) + rttm.Var5(j);

            % IF IT IS LARGER THAN ZERO
            if endTime - startTime > 0

                % GET INDICES FOR START AND END TIME
                startInd = single(find(data.frameTime == startTime));
                endInd   = single(find(data.frameTime == endTime)) - 1;

                % GRAB THIS DATA SUBSET / SEGMENT
                subsetData = data(startInd:endInd,:);

                % Step 4. Grab speech frames that are above speech threshold.
                speechInd = dataVoice(startInd:endInd) > speechCutoff;

                % UPDATE SUBSETDATA
                subsetData = subsetData(speechInd,:);

                % INIT
                thisData_mean   = [];
                thisData_max    = [];
                thisData_min    = [];
                thisData_median = [];
                thisData_std    = [];
                for k = 3:length(names) - 3

                    % GETTING MEAN, MAX, MIN, MEDIAN OF EACH COLUMN
                    % Step 5 (Optional). Calculating segment level statistics
                    if ~contains(names{k},{'F0', 'jitter', 'shimmer'})

                        thisData_mean   = [thisData_mean,   eval(['nanmean(subsetData.',  names{k},')']) ];
                        thisData_max    = [thisData_max,    eval(['nanmax(subsetData.',   names{k},')']) ];
                        thisData_min    = [thisData_min,    eval(['nanmin(subsetData.',   names{k},')']) ];
                        thisData_median = [thisData_median, eval(['nanmedian(subsetData.',names{k},')']) ];
                        thisData_std    = [thisData_std,    eval(['nanstd(subsetData.',   names{k},')']) ];

                        % FOR F0, SHIMMER, and JITTER WE ONLY WANT TO
                        % CALCULATE WITH THE NON-ZEROS
                    else

                        % GETTING PITCH/JITTER/SHIMMER FRAMES WITH
                        % ACTUAL DATA
                        keepInd = eval(['subsetData.',names{k},' ~= 0']);

                        % SUBSTITUTING ZEROS FOR NANs
                        % thowing out zero pitch measures.
                        eval(['subsetData.',names{k},'(~keepInd) = nan * ones(sum(~keepInd),1);']);

                        % NOW DOING THE CALCULATION
                        thisData_mean   = [thisData_mean,   eval(['nanmean(subsetData.',  names{k},')']) ];
                        thisData_max    = [thisData_max,    eval(['nanmax(subsetData.',   names{k},')']) ];
                        thisData_min    = [thisData_min,    eval(['nanmin(subsetData.',   names{k},')']) ];
                        thisData_median = [thisData_median, eval(['nanmedian(subsetData.',names{k},')']) ];
                        thisData_std    = [thisData_std,    eval(['nanstd(subsetData.',   names{k},')']) ];

                    end

                end

                % APPEND SUBJECT DATA
                % the zero pitch measures have been thrown out
                subsetData_spkr0 = [subsetData_spkr0; subsetData];           

            else
                thisData_mean(1:48)     = nan * ones(48,1);
                thisData_max(1:48)      = nan * ones(48,1);
                thisData_min(1:48)      = nan * ones(48,1);
                thisData_median(1:48)   = nan * ones(48,1);
                thisData_std(1:48)      = nan * ones(48,1);
            end
            
        % ADD THE MEANS OF THE SEGMENT TO THE STRUCT
        % -------------------------------------------
        % this is for keeping track of segment level statistics. But we aren't saving it.
        
        % DATA.audio_mean(count,1:48)         = thisData_mean;
        % DATA.audio_max(count,1:48)          = thisData_max;
        % DATA.audio_min(count,1:48)          = thisData_min;
        % DATA.audio_median(count,1:48)       = thisData_median;
        % DATA.audio_std(count,1:48)          = thisData_std;
        
        % count = count+1;


        end

        % NOW MEAN PER SUBJECT/TESTER
        % ----------------------------
        % Step 6. Let's calculate the speaker statistics across all their frames in the file.
        % this is the mean, max, min, median and std.

        thisData_spkr0_mean          = [];
        thisData_spkr0_max           = []; 
        thisData_spkr0_min           = []; 
        thisData_spkr0_median        = []; 
        thisData_spkr0_std           = []; 

        for k = 3:length(names) - 3

            % GETTING MEAN OF EACH COLUMN
            thisData_spkr0_mean = [thisData_spkr0_mean, eval(['nanmean(subsetData_spkr0.',  names{k},')']) ];

            % GETTING MAX OF EACH COLUMN
            thisData_spkr0_max = [thisData_spkr0_max, eval(['nanmax(subsetData_spkr0.',  names{k},')']) ];

            % GETTING MIN OF EACH COLUMN
            thisData_spkr0_min = [thisData_spkr0_min, eval(['nanmin(subsetData_spkr0.',  names{k},')']) ];

            % GETTING MEDIAN OF EACH COLUMN
            thisData_spkr0_median = [thisData_spkr0_median, eval(['nanmedian(subsetData_spkr0.',  names{k},')']) ];

            % GETTING STD OF EACH COLUMN
            thisData_spkr0_std = [thisData_spkr0_std, eval(['nanstd(subsetData_spkr0.',  names{k},')']) ];

        end
        
        % append the subject id and data to the DATA table.
        DATA.name(i,1) = {name};
        
        % append the data
        for NM = {'mean', 'max' , 'min', 'median', 'std'}
            eval(['DATA.audio_spkr0_',char(NM),'{i}  = thisData_spkr0_',char(NM),';']);
        end
            
    end
    
    %% SAVE DATA
    save('DATA.mat','DATA')

    % expand data
    X_data = reshape([DATA.audio_spkr0_mean{:} ...
                      DATA.audio_spkr0_max{:} ...
                      DATA.audio_spkr0_min{:} ...
                      DATA.audio_spkr0_median{:} ...
                      DATA.audio_spkr0_std{:}],44*5,[])';

    % save subjID as num and data
    X_data = [cellfun(@str2num,DATA.name) X_data];
    
    % get new headers
    headers = names(3:end-3);
    X_headers = [];
    for NM = {'mean', 'max' , 'min', 'median', 'std'};
        for kk = 1:length(headers)
            X_headers{end+1} = [headers{kk},'_',char(NM{1})];
        end
    end

    % convert to table
    audio_table = array2table(X_data,'VariableNames',['subjID' X_headers]);

    % save to csv
    writetable(audio_table,'audioFeatures.csv');
    
end

% eof