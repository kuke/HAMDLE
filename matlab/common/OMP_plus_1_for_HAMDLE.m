function [x,it] = OMP_plus_1_for_Legendre(A,b,eps,K1)
% b = b/norm(b); 
% for k=1:size(A,2);
%     A(:,k)=A(:,k)/norm(A(:,k));
% end

ColNorms = zeros(size(A,2),1); 
for k=1:length(ColNorms);
    ColNorms(k)=norm(A(:,k));
end

[M,N] = size(A);

x = zeros(N,1);
r = b;

%A is already normalized
%Anorms = sqrt(sum(A.^2,1)); %Pre compute norms of columns for speed

T = [];
it=0;
z=1;

delta = 1;
% while ((it<K) && (z>0));
while ( ((abs(sum(x)-1) > eps)||delta > 1e-5)  && (it < K1) )
    e = A'*r;  %It is actually faster than e(~T) = A(:,~T)'...;
    e(T)=0;
    e(find(e <= 0))=0;
    e = e ./ ColNorms;
    
    [z,j0]=max(abs(e));
    T = [T j0];
    
    Ap = A(:,T);
    x(T)=lsqnonneg(Ap,b);
    r = b - Ap*x(T);

    it=it+1;      
    delta = norm(r);
end

%fprintf('residual = %f\n', norm(r));


