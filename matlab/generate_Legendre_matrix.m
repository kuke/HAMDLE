%% generate Legendre matrix
function Legendre_mat = generate_Legendre_matrix(order, X)
    m = size(X, 2);
    Legendre_mat = zeros(order+1, m);
    % order = 0
    Legendre_mat(1, :) = ones(1, m);
    % order = 1
    Legendre_mat(2, :) = X;
    % order >= 2
    for n = 3:order+1
        Legendre_mat(n,:) = (2*n-3)*X.*Legendre_mat(n-1,:)/(n-1) ...
                            -(n-2)*Legendre_mat(n-2,:)/(n-1);
    end

end

