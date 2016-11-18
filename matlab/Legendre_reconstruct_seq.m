%% 
function Legendre_reconstruct_seq(seq, order)
    figure('NumberTitle', 'off', 'Name','Sequence sample')
    plot_sequence(seq);
    L = length(seq);
    X = -1:2/(L-1):1;
    Legendre_mat = generate_Legendre_matrix(order, X);
    Legendre_coeff = Legendre_expansion(Legendre_mat, seq, order);
    coeff_len = size(Legendre_coeff, 2);
    
    A_coeff = Legendre_coeff(1:coeff_len/4);
    G_coeff = Legendre_coeff(coeff_len/4+1:coeff_len/2);
    C_coeff = Legendre_coeff(coeff_len/2+1:coeff_len*3/4);
    T_coeff = Legendre_coeff(coeff_len*3/4+1:end);
    
    A_subseq = A_coeff*Legendre_mat;
    G_subseq = G_coeff*Legendre_mat;
    C_subseq = C_coeff*Legendre_mat;
    T_subseq = T_coeff*Legendre_mat;
    
    figure('NumberTitle', 'off', 'Name','Sequence reconstructed')
    subplot(411);  stem(X, A_subseq); axis([-1 1 -0.5 1.5]); ylabel('A');
    set(gca,'FontSize',18);
    subplot(412);  stem(X, G_subseq); axis([-1 1 -0.5 1.5]); ylabel('G');
    set(gca,'FontSize',18);
    subplot(413);  stem(X, C_subseq); axis([-1 1 -0.5 1.5]); ylabel('C'); 
    set(gca,'FontSize',18);
    subplot(414);  stem(X, T_subseq); axis([-1 1 -0.5 1.5]); ylabel('T');
    set(gca,'FontSize',18);
    
    figure('NumberTitle', 'off', 'Name','Reconstructed vs. Original')
    plot_sequence(seq);
   
    subplot(411);  plot(X, A_subseq,'r','lineWidth',2.0); axis([-1 1 -0.5 1.5]); ylabel('A');
    set(gca,'FontSize',18);
    subplot(412);  plot(X, G_subseq,'r','lineWidth',2.0); axis([-1 1 -0.5 1.5]); ylabel('G');
    set(gca,'FontSize',18);
    subplot(413);  plot(X, C_subseq,'r','lineWidth',2.0); axis([-1 1 -0.5 1.5]); ylabel('C');
    set(gca,'FontSize',18);
    subplot(414);  plot(X, T_subseq,'r','lineWidth',2.0); axis([-1 1 -0.5 1.5]); ylabel('T');
    set(gca,'FontSize',18);
end