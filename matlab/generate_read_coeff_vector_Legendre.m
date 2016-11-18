%% 
function generate_read_coeff_vector_Legendre(READ_seq, Legendre_mat, order)

    L = size(Legendre_mat, 2);

    if order+1 > size(Legendre_mat, 1)
        X = -1:1/(L-1):1;
        Legendre_mat = generate_Legendre_matrix(order, X);
    end
 
    reads = READ_seq;
    
    n = size(reads, 1);

    Legendre_coeff_vector_reads = zeros(n, 4*(order+1));
   
    disp('Generating Legendre coeff vector for reads');

    for i = 1:n
        seq = reads(i, 1).Sequence;
        Legendre_coeff_vector_reads(i,:) = Legendre_expansion(Legendre_mat, seq(1:L), order);
    end
    
    filename = 'Legendre_coeff_vector_reads.mat';
    save(filename,'Legendre_coeff_vector_reads');
end