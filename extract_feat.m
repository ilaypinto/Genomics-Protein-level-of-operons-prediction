function features = extract_feat(data)
    % preallocate memory
    features = zeros(size(data,1),51);
    couples = ["GG","GC","GA","GT","CC","CG","CA","CT","AA","AG","AC","AT","TT",...
    "TG","TC","TA"];             % All couple nucleotides options

    for i = 1:size(data,1)
    
        % Relevant data
        NT = data{i,2}{:}(27:50);
        AA = nt2aa(data{i,2}{:}(27:50));
        
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
        end
        features(i,28) = length(strfind(NT,'TATA')); % TATA is in Eukaryotes so maybe no...
        features(i,29) = length(strfind(NT,'CAT'));  % Saw in wikipedia, idk...
        features(i,30) = length(strfind(NT,'A'))+length(strfind(NT,'G')); % No. of Purines
        features(i,31) = length(strfind(NT,'C'))+length(strfind(NT,'T')); % No. of Pyrimidines
    
        % Task 1 - AAA,TTT, GCA
        features(i,32) = length(strfind(NT,'AAA'));
        features(i,33) = length(strfind(NT,'TTT'));
        features(i,34) = length(strfind(NT,'GCA'));
    
        features(i,35) = rnafold(NT);
        for j = 1:length(couples)
            features(i,35+j) = length(strfind(NT,couples(1,j)));
        end
        
        % next feature is no 52 (35+16 couples)
        
    end
end
