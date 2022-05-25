%% Main Script - Challenge %%
%#ok<*UNRCH> 
clear all; clc; close all;

%% Flags
% add flags as you please for efficiant work
create_bins = 1;         % 1 - Convert to bin based features, 0 - dont convert into bin based features
feat_selection = 1;      % 1 - apply feature selection process, 0 - dont aplly feature selection
val_flag = 1;            % 1 - Prediction and analysis on Validation set flag
test_flag = 1;           % 1 - Prediction and analysis on Test set flag
finito_flag = 1;         % 1 - Prediction on Unknown Data flag
use_all_feat = 1;        % 1 - use all features, 0 - use only features from xlsx
extract_feat_flag = 1;        % 1 - extract features from raw data, 0 - load features mat

% Useful Parameters
    couples = ["GG","GC","GA","GT","CC","CG","CA","CT","AA","AG","AC","AT","TT",...
        "TG","TC","TA"];  % Maybe use this?

% Load Data
all_data  = readtable('Variants_sequence.xlsx', 'VariableNamingRule', 'preserve');   % Load all sequences

%% Features Creation
if extract_feat_flag 
    features = extract_feat(all_data);
    save('features.mat', "features");
else 
    features = load('features.mat');
    features = features.features;
end


%% Convert to bin based features (mean of features from sequences in the bin)
if use_all_feat
    [X_known_our, targets_1, X_unknown_our] =  bin_feat(features);  
    [X_known_xlsx, X_unknown_xlsx, targets] = bin_feat_from_files();
    
    X_known   = cat(2, X_known_xlsx, X_known_our);
    X_unknown = cat(2, X_unknown_xlsx, X_unknown_our);

% quick sanity check
    if targets_1 ~= targets
        disp('oh no this is bad! there is a miss match between the two target vectors...');
    end
else
    [X_known, X_unknown, targets] = bin_feat_from_files(); 
end


%% Correlation tests and feature rejection
feat_feat_corr = abs(corr(X_known,type = 'Spearman'));
feat_target_corr = abs(corr(X_known,targets,type = 'Spearman'));
% insert here the corr analysis function!

%% Feature selection
X_known = rand(300,30); targets = rand(300,1);
if feat_selection  
    idx = randperm(length(targets), length(targets)*0.5); % use 50% from the data to select features (prevent overfiting)
    % remove more features via SFS or filter methods.
    options = statset('Display', 'iter', 'UseParallel', true);  % UseParallel to speed up the computations and Display so we can see the progress
    fun = @(Xtrain,Ytrain,Xtest,Ytest) sfs_corr(Xtrain,Ytrain,Xtest,Ytest);  % sfs under correlation of linear regression predictions with targets
    [idx_sfs, history_sfs] = sequentialfs(fun, X_known(idx,:), targets(idx), 'options', options); 
   
    save('idx_sfs.mat', 'idx_sfs');
else
    idx_sfs = load('idx_sfs.mat');
    idx_sfs = idx_sfs.idx_sfs;
end

X_known(:,~idx_sfs) = [];
X_unknown(:,~idx_sfs) = [];

% figure;
% gplotmatrix(X_known, [], targets);  % gplot - look at the features!


%% Prediction and analysis using a 10-fold cross validation
C_2 = cvpartition(length(targets), 'KFold', 10); % cvpartition object - creating a 10-fold partition

% train a model - important to add catagorical features indices (if there are any) when training the model 
linear_model = fitrlinear(X_known,  targets, "CVPartition", C_2, "Learner", "leastsquares"); % Simple linear regression
svm_model    = fitrsvm(X_known,     targets, "CVPartition", C_2, "KernelFunction","linear"); % SVM regression
neural_model = fitrnet(X_known,     targets, 'CVPartition', C_2); % Neural Network regression
tree_model   = fitrensemble(X_known,targets, "CVPartition", C_2);  % regression ensemble Tree

% predict targets using the partitioned models
Y_linear = kfoldPredict(linear_model);
Y_svm    = kfoldPredict(svm_model);
Y_NN     = kfoldPredict(neural_model);
Y_tree   = kfoldPredict(tree_model);

% compute correlation between prediction and targets - dont use abs cause
% we are also interest in a positive correlation!
correlations = [
    corr(targets, Y_linear, type = 'Spearman')
    corr(targets, Y_svm,    type = 'Spearman')
    corr(targets, Y_NN,     type = 'Spearman')
    corr(targets, Y_tree,   type = 'Spearman')];

% visualization
figure('Name', 'regression plot')
plot(targets); hold all; 
plot(Y_linear); % linear model
plot(Y_svm); % svm model
plot(Y_NN); % NN model
plot(Y_tree); % tree model
title('model comparison'); legend({'targets', 'linear predictions', 'svm predictions', 'NN predictions', 'tree predictions' });

%% train a model on the entire data
% train a model - important to add catagorical features indices (if there are any) when training the model 
final_linear_model = fitrlinear(X_known,   targets, "Learner", "leastsquares"); % Simple linear regression
final_svm_model    = fitrsvm(X_known,      targets); % SVM regression
final_NN_model     = fitrnet(X_known,      targets); % Neural Network regression
final_tree_model   = fitrensemble(X_known, targets); % regression Tree

%% Prediction on Unknown Data
if finito_flag
    predict_test_linear = predict(final_linear_model, X_unknown);
    predict_test_svm    = predict(final_svm_model,    X_unknown);
    predict_test_NN     = predict(final_NN_model,     X_unknown);
    predict_test_tree   = predict(final_tree_model,   X_unknown);
end


%% lets see how close we are on the unknown (thank you data leaks :))
% define the paths for the xlsx files and load the data from them
% known_path = 'known_data_set.xlsx';
% unknown_path = 'unknown_data_set.xlsx';
% 
% known = sortrows(readtable(known_path, "VariableNamingRule", "preserve"), 'Bin index') ;
% unknown = sortrows(readtable(unknown_path, "VariableNamingRule", "preserve"), 'Bin index');
% 
% % extract the bin numbers
% known_bin_idx = known(:,1).Variables;     
% unknown_bin_idx = unknown(:,1).Variables;
% 
% [sorted, I] = sort([known_bin_idx ;unknown_bin_idx]);
% t = [targets ;predict_test_svm];
% 
% plot(sorted, t(I));

% C:\Users\tomer\Documents\ViennaRNA-2.5.0
