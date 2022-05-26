function features = extract_feat(data)
    % preallocate memory
    features = zeros(size(data,1),95);
    couples = ["GG","GC","GA","GT","CC","CG","CA","CT","AA","AG","AC","AT","TT",...
    "TG","TC","TA"];             % All couple nucleotides options
    f = waitbar(0, ['pls wait! ' num2str(0) ' out of ' num2str(size(data,1))]);
    for i = 1:size(data,1)
        waitbar(i/size(data,1), f, ['pls wait! ' num2str(i) ' out of ' num2str(size(data,1))]); % update the wait bar

        % Relevant data
        NT = data{i,2}{:};
        AA = nt2aa(NT);
        
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
        features(i,32) = length(strfind(NT,'TAATG'));                     % Sequence from paper
        features(i,33) = length(strfind(NT,'TGATG'));                     % Sequence from paper
        features(i,34) = length(strfind(NT,'TAGTG'));                     % Sequence from paper
        features(i,35) = length(strfind(NT,'TAGA'));                      % Sequence from paper
    
        % Task 1 - AAA,TTT, GCA
        features(i,36) = length(strfind(NT,'AAA'));
        features(i,37) = length(strfind(NT,'TTT'));
        features(i,38) = length(strfind(NT,'GCA'));
    
        [~,features(i,39)] = rnafold(NT);
        
        for j = 1:length(couples)
            features(i,39+j) = length(strfind(NT,couples(1,j)));
        end

        win_size = 30;
        for j = 1:(length(NT) - win_size)
%             system(['echo "' NT(j:j+win_size) '" >> C:\Users\tomer\Documents\Fold\sequences']);
            [~,features(i,55 + j)] = rnafold(NT(j:j+win_size));
        end
%         [~,features(i,56:end)] = system(['C:\Users\tomer\Documents\Fold\RNAfold < C:\Users\tomer\Documents\Fold\sequences']);
%         delete 'C:\Users\tomer\Documents\Fold\sequences'
    end
    delete(f)
end
