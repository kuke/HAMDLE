function [ C, w ] = LBG(T, C, w)
% T is the data to be clustered, where rows represent datapoints and
% columns the dimension of the data. C is the matrix with the centroids of
% each cluster, here row represents the cluster and the columns are the
% dimensions of the cluster centroid. W is the corresponding weights of
% each cluster.

eps = 0.000005;
[M, k] = size(T);
N = length(w);

% Squared-error distortion measure.
dist = 1/(M*k)*sum(pdist2(C, T, 'euclidean', 'smallest', 1).^2);
dist_old = 0;

% Find index of the largest cluster.
[~, I] = max(w);

% Split it
C = [C; zeros(1, k)];
C(I, :) = (1+eps)*C(I, :);
C(N+1, :) = (1-eps)*C(I, :);
N = N + 1;

iter = 0;
delta = 1;

while delta > eps
    % D is distance and Q is its index in C, returns only the smallest
    % euclidean distances.
    [D, Q] = pdist2(C, T, 'euclidean', 'Smallest', 1);
    
    w = zeros(1, N);
    C = zeros(N, k);
    
    % Update the weights and recompute centroid of each cluster
    for i = 1:length(Q)
        w(Q(i)) = w(Q(i)) + 1;
        C(Q(i),:) = C(Q(i),:) + T(i,:);
    end
    for i = 1:length(w)
        if w(i) == 0
            C(i,:) = 0;
        else
            C(i,:) = C(i,:)/w(i);
        end
    end
    
    w = w/sum(w);
    iter = iter + 1;
    dist_old = dist;
    dist = 1/(M*k)*sum(D.^2);
    delta = (dist - dist_old)/dist_old;
end
end