function [known_feat, unknown_feat, known_labels] = bin_feat_from_files()

% define the paths for the xlsx files and load the data from them
known_path = 'known_data_set.xlsx';
unknown_path = 'unknown_data_set.xlsx';

known = sortrows(readtable(known_path, "VariableNamingRule", "preserve"), 'Bin index') ;
unknown = sortrows(readtable(unknown_path, "VariableNamingRule", "preserve"), 'Bin index');

% extract labels for known bins 
known_labels = known(:,28).Variables;

% extract features
known_feat = known(:,29:end).Variables;
unknown_feat = unknown(:,28:end).Variables;

end