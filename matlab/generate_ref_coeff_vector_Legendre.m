function generate_ref_coeff_vector_Legendre(REF_seq, Legendre_mat, order)
    L = size(Legendre_mat, 2);
    Lp = 1;
    if order+1 > size(Legendre_mat, 1)
        X = -1:1/(L-1):1;
        Legendre_mat = generate_Legendre_matrix(order, X);
    end

    reference_seqs = REF_seq;
    
    n = size(reference_seqs, 1);
    all_species_Legendre_coeff_vector = cell(n, 1);

    for i = 1:n
        N = length(reference_seqs(i,1).Sequence)-L+1;
        seq = reference_seqs(i,1).Sequence;
        species_Legendre_coeff_vector = zeros(N, 4*(order+1));

        for j = 1:N
               species_Legendre_coeff_vector(j,:) = Legendre_expansion(Legendre_mat, seq(j:(j+L-1)), order);
        end
        all_species_Legendre_coeff_vector{i,1} = species_Legendre_coeff_vector;
    end 

    filename = strcat('Legendre_coeff_vector_allref.mat');
    save(filename,'all_species_Legendre_coeff_vector');
end 