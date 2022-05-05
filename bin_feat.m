function [known_bin_idx_sorted, known_bin_features, known_labels, unknown_bin_idx_sorted,...
          unknown_bin_features] =  bin_feat(features_mat)

% define the paths for the xlsx files and load the data from them
known_path = 'known_data_set.xlsx';
unknown_path = 'unknown_data_set.xlsx';

known = readtable(known_path, "VariableNamingRule", "preserve");
unknown = readtable(unknown_path, "VariableNamingRule", "preserve");

% extract the bin numbers
known_bin_idx = known(:,1).Variables;     
unknown_bin_idx = unknown(:,1).Variables;

% extract the indices of the sequences for each bin
known_seq_idx = known(:,2:27).Variables;
unknown_seq_idx = unknown(:,2:27).Variables;

% extract labels for known bins 
known_labels = known(:,28).Variables;

% initialize an empty matrix for the bin features
bin_features = zeros(max([known_bin_idx;unknown_bin_idx]), size(features_mat, 2)+10);

% compute the bin features and allocate it into its correct idx in the bin
% features matrix
% start with the known bins
for i = 1:size(known_bin_idx)
    idx = known_bin_idx(i);                      % bin index
    seq_indices = known_seq_idx(i,:);            % sequences indices in the bin
    seq_indices(isnan(seq_indices)) = [];
    curr_features = features_mat(seq_indices,:); % features of that sequences
    curr_bin_Feat = mean(curr_features);         % mean of each column (feature)
    bin_features(idx,:) = [curr_bin_Feat,table2array(known(i,29:38))];% allocate the features into the features matrix
end

% repeat the same loop only with the unknown bins
for i = 1:size(unknown_bin_idx)
    idx = unknown_bin_idx(i);                      % bin index
    seq_indices = unknown_seq_idx(i,:);            % sequences indices in the bin
    seq_indices(isnan(seq_indices)) = [];
    curr_features = features_mat(seq_indices,:); % features of that sequences
    curr_bin_Feat = mean(curr_features);         % mean of each column (feature)
    bin_features(idx,:) = [curr_bin_Feat,table2array(unknown(i,28:37))];         % allocate the features into the features matrix
end

% sort the bin indices to keep alignment between bins nambers and indices
% in the 'bin_features' matrix
[known_bin_idx_sorted, I] = sort(known_bin_idx);
unknown_bin_idx_sorted = sort(unknown_bin_idx);


% extract each set features (known\unknown)
known_bin_features = bin_features(known_bin_idx_sorted,:);
unknown_bin_features = bin_features(unknown_bin_idx_sorted,:);

known_labels = known_labels(I); % rearange the labels of known to match the feature matrix

end

