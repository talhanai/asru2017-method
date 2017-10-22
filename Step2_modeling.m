%% EVALUTING FEATURES USING REGULARIZED LOG REG
%
% Author: Tuka Alhanai, CSAIL MIT 2017 Copyright

clear all

% you will need to install the glmnet package
addpath(genpath('../tools/glmnet_matlab'))

% tunring off some annoying warnings.
warning('off','stats:glmfit:PerfectSeparation')
warning('off','stats:glmfit:IterationLimit')
warning('off','stats:glmfit:IllConditioned')

 
% LOAD Features

% outcomes is 2 columns, col1 is subjectID, col2 is outcome 1/0.
outcomes  = readtable('outcomes.csv');

% features have first column with subject ID, rest of the columns are features.
audioData = readtable('audioFeatures.csv');

% load acoustic features
X_data = zscore(table2array(audioData));

% grabbing names of audio features
X_names = audioData.Properties.VariableNames;


% throwout first column which is subject ID.
X_data(:,1) = [];
X_names(1) = [];

% grab subject labels
labels = outcomes.label;

% grab subject IDs
subjID = outcomes.subjID;


%% MODELING

% start parallel process
parpool(9)

% init vars
clear cal cal_train
allCorr = []; topCorr = [];
allCorr_dev = []; topCorr_dev = [];
pT = []; allp = [];
glm = []; glmnet_prob = [];
count = 1;
corrT = 0;

% defining p-val threshold, 
for pT = [0.05 0.01 0.001]

    % looping through multiple alphas (L1/L2 penalties)
    for alpha = [0:0.01:1]
    
    % printing out vals
    [pT corrT alpha]
    
        % for reproducibility
        rng('default')

        % subjID = labels(:,1);
        allCorr = [];

        % defining glmnet parameter
        options.alpha = alpha;

        % PERFORM LEAVE ONE OUT CROSS VALIDATION
        for i = 1:size(labels,1)
            
            i

            % GET TEST/TRAIN INDICES
            testInd  = ismember(subjID,subjID(i));
            trainInd = ~testInd;

            % GET TRAIN/TEST SET
            X_test  = X_data(testInd,:);
            X_train = X_data(trainInd,:);
            
            % GET LABELS
            labels_true_test  = labels(testInd);
            labels_true_train = labels(trainInd);
            
            % GET FEATURES THAT ARE CORRELATED WITH LABELS ABOVE THRESHOLD
            [topCorr,p] = corr(X_train,labels_true_train);
            indp = (p < pT);
            allCorr = [allCorr; topCorr'];
            allp    = [allp; p];
            
            % THROW OUT FEATURES ABOVE p THRESHOLD
            X_train(:,~indp) = [];
            X_test(:,~indp)  = [];
        
            % GLM
            % -----
            model = fitglm(X_train,labels_true_train,'Distribution','binomial','link','logit');
            
            glm(i).prob_test   = predict(model, X_test);
            glm(i).label_test  = labels_true_test;
            glm(i).prob_train  = predict(model, X_train);
            glm(i).label_train = labels_true_train;
            glm(i).indp        = indp;
            
            % GLMNET MODEL
            % -------------
            CVfit = cvglmnet(X_train, labels_true_train, 'binomial',options,'deviance',size(X_train,1),[],true);
            
            glmnet_prob(i).test  = cvglmnetPredict(CVfit,X_test,[],'response');
            glmnet_prob(i).train = cvglmnetPredict(CVfit,X_train,[],'response');
            glmnet_prob(i).labels_test  = labels_true_test; 
            glmnet_prob(i).labels_train = labels_true_train; 

        end

        % GET PERFORMANCE
        
        % glm performance
        [~,~,~,AUC]              = perfcurve([glm(:).label_test] , [glm(:).prob_test]  , 1);
        cal.test(count)          = AUC;
        [~,~,~,AUC]              = perfcurve( reshape([glm(:).label_train],[],1) , reshape([glm(:).prob_train],[],1)  , 1);
        cal.train(count)         = AUC;
    
        % variable data
        cal.corrInd{count}       = allCorr;
        cal.corrT(count)         = corrT;
        cal.pT(count)            = pT;
        cal.indp{count}          = [glm(:).indp];
        cal.alpha(count)         = alpha;
                
        % glmnet performance
        [~,~,~,AUC]              = perfcurve([glmnet_prob(:).labels_test] , [glmnet_prob(:).test]  , 1);
        cal.test_glmnet(count)   = AUC;
        [~,~,~,AUC]              = perfcurve( reshape([glmnet_prob(:).labels_train],[],1) , reshape([glmnet_prob(:).train],[],1)  , 1);
        cal.train_glmnet(count)  = AUC;

        count                    = count+1 ;

    end

end