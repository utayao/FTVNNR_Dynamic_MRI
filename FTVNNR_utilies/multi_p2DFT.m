function A = multi_p2DFT(mask)
%
%   multi_p2DFT
%
[m, n, T] = size(mask); 
A_single = cell(1,T);
for i = 1 : T
   A_single{i} = p2DFT(mask(:,:,i), [m n], 1, 2); 
end
 %  A = A_operator(@(X) multi_A(A_single, X, T, [m, n]), @(X) multi_A(cellfun(@(X) X',A_single), X, T, [m, n])); 
 A = A_operator(@(X) multi_A(A_single, X, T, [m, n]), @(X) multi_B(A_single, X, T, [m, n])); 
end

function R = multi_A(A_single, X, T, X_single_size)
    R = zeros(size(X));
    N = prod(X_single_size);
    for t = 1:T
        idx = (t-1)*N;
        R(idx+1:idx+N) = A_single{t}*reshape(X(idx+1:idx+N), X_single_size);
    end
end

function R = multi_B(A_single, X, T, X_single_size)
    R = zeros(size(X));
    N = prod(X_single_size);
    for t = 1:T
        idx = (t-1)*N;
        R(idx+1:idx+N) = A_single{t}'*reshape(X(idx+1:idx+N), X_single_size);
    end
end
