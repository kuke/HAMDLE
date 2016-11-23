%% compute between-class max cos similarity
function max_sim = compute_between_class_max_cos_similarity(A, B)
    [m_A, n_A] = size(A);
    [m_B, n_B] = size(B);
    norm_Ai = sqrt(sum(A.*A, 1));
    for i = 1:n_A
        A(:, i) = A(:, i)/norm_Ai(i);
    end
    norm_Bi = sqrt(sum(B.*B, 1));
    for i = 1:n_B
        B(:, i) = B(:, i)/norm_Bi(i);
    end
    
    max_sim = max(max(A'*B));
  