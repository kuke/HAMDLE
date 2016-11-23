function [result_Legendre elapsedtime_READ elapsedtime_REF elapsedtime_algo it] = Legendre_implementation(order,Legendre_mat, READ_seq,REF_seq, NoOfSpecies, seq2species)

    I=1000; % Maximum allowable iteration in OMP^{+,1}    
    Nu=0.00001; % allowable tolerance in l1 norm computation for OMP^{+,1} 
    %
    tic;
    tstart_READ = tic;
    generate_read_coeff_vector_Legendre(READ_seq, Legendre_mat, order);
    elapsedtime_READ = toc(tstart_READ);
    
    fprintf('Elapsed time to generate Legendre coeff vector for reads: elapsedtime_READ = %f\n', elapsedtime_READ);
 
    tic;
    tstart_REF = tic;
    generate_ref_coeff_vector_Legendre(REF_seq, Legendre_mat, order);
    elapsedtime_REF = toc(tstart_REF);
    fprintf('Elapsed time to generate Legendre coeff vector for reference: elapsedtime_REF = %f\n', elapsedtime_REF);
    
    allref_vector_compile_Legendre('Legendre_coeff_vector_allref.mat', seq2species);
    
    %% SEK (using OMP^{+,1}) execution

    % Loading the relevant mat files
    load('Legendre_coeff_vector_reads.mat'); load('Legendre_coeff_vector_allref_trans_Legendre.mat'); 
    X = all_species_Legendre_coeff_trans;   
    Mu = mean(Legendre_coeff_vector_reads,1)';

    tstart_algo=tic;
    [gamma it]= OMP_plus_1_for_Legendre(X,Mu,Nu,I); gamma = gamma / sum(gamma);
    fprintf('Iteration times: it = %f\n', it);

    result_Legendre = zeros(1,NoOfSpecies);

    for i =1:length(gamma)
        if gamma(i) ~= 0, result_Legendre(fragment2species(i)) = result_Legendre(fragment2species(i)) + gamma(i); end
    end
    elapsedtime_algo = toc(tstart_algo);

    fprintf('Elapsed time to execute OMP^{+,1} algorithm: elapsedtime_algo = %f\n', elapsedtime_algo); 
end