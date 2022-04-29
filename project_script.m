%% Main Script - Challenge %%

% This is the main script for the project 
clear all; clc; close all;

%% Flags
% add flags as you please for efficiant work
feat_creation_flag = 1;  % Features Creation flag
convert_flag = 1;        % Convert to bin based features flag
feat_selection_flag = 1; % Feature selection flag
eval_flag = 1;           % Prediction and analysis on Validation set flag
test_flag = 1;           % Prediction and analysis on Test set flag
finito_flag = 1;         % Prediction on Unknown Data flag


% Useful Parameters
    couples = ["GG","GC","GA","GT","CC","CG","CA","CT","AA","AG","AC","AT","TT",...
        "TG","TC","TA"];  % Maybe use this?

%% Features Creation
if feat_creation_flag == 1
    % Load Data
    all_data  = readtable('Variants_sequence.xlsx');   % Load all sequences
    known_data = readtable("known_data_set.xlsx");     % Load known data
    unknown_data = readtable("unknown_data_set.xlsx"); % Load unknown data
    
    features = zeros(size(all_data,1),50);

    for i = 1:size(all_data,1)
    
        % Relevant data
        NT = all_data{i,2}{:}(27:50);
        AA = nt2aa(all_data{i,2}{:}(27:50));
        
        % Task 3 - more features!
    
        % Amino Acid based features
        features(i,1) = length(strfind(AA,'A')); % Alanine
        features(i,2) = length(strfind(AA,'R')); % Arganine
        features(i,3) = length(strfind(AA,'N')); % Asparagine
        features(i,4) = length(strfind(AA,'D')); % Asparatate
        features(i,5) = length(strfind(AA,'C')); % Cysteine
        features(i,6) = length(strfind(AA,'Q')); % Glutamine
        features(i,7) = length(strfind(AA,'E')); % Glutamate
        features(i,8) = length(strfind(AA,'G')); % Glycine
        features(i,9) = length(strfind(AA,'H')); % Histidine
        features(i,10) = length(strfind(AA,'I')); % Isoleucine
        features(i,11) = length(strfind(AA,'L')); % Leucine
        features(i,12) = length(strfind(AA,'K')); % Lysine
        features(i,13) = length(strfind(AA,'M')); % Methionine
        features(i,14) = length(strfind(AA,'F')); % Phenylalanine
        features(i,15) = length(strfind(AA,'P')); % Proline
        features(i,16) = length(strfind(AA,'S')); % Serine
        features(i,17) = length(strfind(AA,'T')); % Threonine
        features(i,18) = length(strfind(AA,'W')); % Tryptophan
        features(i,19) = length(strfind(AA,'Y')); % Tyrosine
        features(i,20) = length(strfind(AA,'V')); % Valine
        features(i,21) = length(strfind(AA,'*')); % Stop codon
        features(i,22) = features(i,2)+features(i,9)+features(i,12)+features(i,4)+...
          features(i,7);                        % Electricaly Charged side chains
        features(i,23) = features(i,16)+features(i,17)+features(i,3)+features(i,7)+...
          features(i,6);                        % Polar Uncharged side chains
        features(i,24) = features(i,1)+features(i,20)+features(i,10)+features(i,11)+...
          features(i,13)+features(i,14)+features(i,19)+features(i,18);
                                                % Hydrophobic side chains
        features(i,25) = features(i,5)+features(i,8)+features(i,15);
                                                % Special cases
        % DNA based features
        [pal_pos,pal_length] = palindromes(NT);
        if ~isempty(pal_pos)
            features(i,26) = length(pal_pos);       % Number of palindromes in sequence 
            features(i,27) = max(pal_length);       % Longest palindrome in sequence
        else
            features(i,26) = 0;       % Number of palindromes in sequence
            features(i,27) = 0;       % Longest palindrome in sequence
            % features 26-27 are highly corroloated, consider removing one.
        end
        features(i,28) = length(strfind(NT,'TATA')); % TATA is in Eukaryotes so maybe no...
        features(i,29) = length(strfind(NT,'CAT'));  % Saw in wikipedia, idk...
        features(i,30) = length(strfind(NT,'A'))+length(strfind(NT,'G')); % No. of Purines
        features(i,31) = length(strfind(NT,'C'))+length(strfind(NT,'T')); % No. of Pyrimidines
            % features 30-31 are highly corroloated, consider removing one.
    
    
        % Task 1 - AAA,TTT, GCA
        features(i,32) = length(strfind(NT,'AAA'));
        features(i,33) = length(strfind(NT,'TTT'));
        features(i,34) = length(strfind(NT,'GCA'));
    
        % Task 2 - Folding energy
    %     features(i,35) = FE window 1
    %     features(i,36) = FE window 2
    %     features(i,37) = FE window 3
    
   
    end
    save('features.mat','features');   % Save created features
else
    load('features.mat');              % Load created features
end
%% Convert to bin based features
if convert_flag == 1
    [known_bin_idx_sorted, X_known, Y, unknown_bin_idx_sorted,...  % Converting created features
    X_unknown] =  bin_feat(features);                              % into bin based features
    
    save('bin_convertion.mat','X_known','Y',...
    'X_unknown','known_bin_idx_sorted','unknown_bin_idx_sorted');  % Save created parameters

else
    load('bin_convertion.mat'); % Load created parameters
end

%% Correlation

heatmap(abs(corr(X_known,type = 'Spearman')))
heatmap(abs(corr(X_knwon,Y,type = 'Spearman')))

%% Feature selection
if feat_selection_flag == 1
    % Before removing features- consider feature-feature and feature label
    % correlation.
    corr_remove = [31,32,33,34];    % features indices to remove due to coorelation analysis
    X_knwon(:,corr_remove) = [];
    X_unknwon(:,corr_remove) = [];
    
    % remove more features via SFS or filter methods.
    options = statset('Display', 'iter', 'UseParallel', true);  % UseParallel to speed up the computations and Display so we can see the progress
    fun = @(Xtrain,Ytrain,Xtest,Ytest)loss(fitlm(Xtrain, Ytrain), Xtest, Ytest);  % sfs under loss of linear regression
    [idx_sfs, history_sfs] = sequentialfs(fun, X_known, Y, 'options', options); 
    
    sfs_remove = ~idx_sfs;
    X_knwon(:,sfs_remove) = [];
    X_unknwon(:,sfs_remove) = [];
    
    save('feat_selection.mat','corr_remove','sfs_remove','X_known','X_unknown');
else
    load('feat_selection.mat');
end
%% Prediction and analysis on Validation set
if eval_flag == 1

    % Create test and training/validation sets
    indices=cvpartition(Y,'holdout',0.2);
    training_indices = training(indices);
    test_indices = ~training_indices;
    features_mid = X_known(training_indices,:);
    labels_mid = Y(training_indices,:); 
    features_test = X_known(test_indices,:);
    labels_test =Y(test_indices,:); 
    
    % Creat training and validation sets
    indices=cvpartition(labels_mid,'holdout',0.5);
    training_indices = training(indices);
    validation_indices = ~training_indices;
    features_training = features_mid(training_indices,:);
    labels_training = labels_mid(training_indices,:); 
    features_validation = features_mid(validation_indices,:);
    labels_validation =labels_mid(validation_indices,:);
    
    % Task 4 - Create model on training set
    linear_model = fitlm(features_training,labels_training);      % Simple linear regression
    neural_model = fitrnet(features_training,labels_training);    % Neural Network regression
    tree_model = fitrensemble(features_training,labels_training); % Trees ensemble regression
    
    % Predict on validation set and on training set itself
    % On validation
    linear_predict_val = predict(linear_model, features_validation);
    neural_predict_val = predict(neural_model, features_validation);
    tree_predict_val = predict(tree_model, features_validation);
    
    % On training 
    linear_predict_train = predict(linear_model, features_training);
    neural_predict_train = predict(neural_model, features_training);
    tree_predict_train = predict(tree_model, features_training);
    
    % Task 5 - Correlation of training vs validation;
    
    Validation_corr = [abs(corr(labels_validation,linear_predict_train,type = 'Spearman')),...
        abs(corr(labels_validation,neural_predict_train,type = 'Spearman')),...
        abs(corr(labels_validation,tree_predict_train,type = 'Spearman'))];
    
    Training_corr = [abs(corr(labels_training,linear_predict_train,type = 'Spearman')),...
        abs(corr(labels_training,neural_predict_train,type = 'Spearman')),...
        abs(corr(labels_training,tree_predict_train,type = 'Spearman'))];
end
%% Prediction and analysis on Test set
if test_flag == 1
    % Create model using the whole training set(training+validation)
    Best_model = fitlm(features_mid,labels_mid); 
    
    % Task 5 - Predict on test set
    predict_test = predict(Best_model, features_test);
end
%% Prediction on Unknown Data
if finito_flag == 1
    % Create model using the whole known set
    Final_model = fitlm(X_knwon,Y); 
    
    % Predict on unknown set
    predict_test = predict(Final_model, X_unknown);

end
