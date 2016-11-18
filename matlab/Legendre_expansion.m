% Legendre expansion
%   @parma
%   Legendre_mat: Legendre matrix
%   seq: input sequence
%   Legendre_mat & seq should be euqal width
function Legendre_coeff = Legendre_expansion(Legendre_mat, seq, order)
    %Legendre_coeff = zeros(1, 4*(order+1));
    lenSeq = size(seq, 2);
    %pos = -1:2/(lenSeq-1):1;

    A_subseq = (seq=='A');
    G_subseq = (seq=='G');
    C_subseq = (seq=='C');
    T_subseq = (seq=='T');

    norm_factor = (1:2:2*order+1)/lenSeq;

    A_coeff = norm_factor.*(A_subseq*Legendre_mat');
    G_coeff = norm_factor.*(G_subseq*Legendre_mat');
    C_coeff = norm_factor.*(C_subseq*Legendre_mat');
    T_coeff = norm_factor.*(T_subseq*Legendre_mat');
    
    Legendre_coeff = [A_coeff, G_coeff, C_coeff, T_coeff];
   
end