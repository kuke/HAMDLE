function allref_vector_compile_Legendre(infile, seq2species)

    % This script is used to compile the original reference vectors into the
    % format which could be processed by compressed sensing script.

    load(infile); % This is the original one.
    %load('GroundTruth.mat'); % This includes the informations about which reference sequence belongs to which species.
    
    nr_column = 0;
    nr_row = 0;
    
    for i = 1:size(all_species_Legendre_coeff_vector,1)
        
        nr_row = size(all_species_Legendre_coeff_vector{i,1},2);
        nr_fragment = size(all_species_Legendre_coeff_vector{i,1},1);
        nr_column = nr_column + nr_fragment;
        
    end
    
    
    all_species_Legendre_coeff_trans = zeros(nr_column, nr_row);
    
    fragment2species = zeros(1,nr_column);
    
    count = 0; nr_column = 0;
    
    for i = 1:size(all_species_Legendre_coeff_vector,1)
        
%         disp(['Right now, No. ' num2str(i) ' species is under processing.']);
        
        nr_fragment = size(all_species_Legendre_coeff_vector{i,1},1);
        nr_column = nr_column + nr_fragment;
        
        for j = 1:nr_fragment
            
            count = count + 1;
            
            all_species_Legendre_coeff_trans(count,:) = all_species_Legendre_coeff_vector{i,1}(j,:);
            
            fragment2species(1, count) = seq2species(i);
            
        end
        
    end
    
    all_species_Legendre_coeff_trans = all_species_Legendre_coeff_trans';
    
    disp(['There are totally ' num2str(nr_column) ' fragments!']);
    
    filename = strcat('Legendre_coeff_vector_allref_trans_Legendre.mat');
    
    save(filename,'all_species_Legendre_coeff_trans','fragment2species');
    
end

    