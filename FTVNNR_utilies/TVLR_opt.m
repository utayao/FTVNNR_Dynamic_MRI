function X = TVLR_opt(A, B, pars,F_gt)
%
%   Arguments
%       A: The downsampling Fourier Transform Operator  
%       B: Multi-channel data : m * n * T * c
%       pars: pars.lambda_1, pars.lambda_2, pars.gamma, pars.max_iter,
%       pars.tol and others
%
%   Output 
%       x: denoised image
%       psnr_vals: psnr vals for every iter
%       energy_vals: energy vals for every iter
%       time_vals: time vals for every iter
%
%   References
%       Jiawen Yao, Zheng Xu, Xiaolei Huang, Junzhou Huang  
%       "An Efficient Algorithm for Dynamic MRI Using Low-Rank and Total Variation Regularizations",
%       In Medical Image Analysis, 44, 14-27, 2018.
%       Jiawen Yao, Zheng Xu, Xiaolei Huang, Junzhou Huang 
%       "Accelerated dynamic MRI reconstruction with total variation and nuclear norm regularization", 
%       In MICCAI 2015.

%% Set Default Parameter
if ~isfield(pars, 'verbose') 
    pars.verbose = 0;
end
if ~isfield(pars, 'debug_output')
    pars.debug_output = 0;
end
if ~isfield(pars, 'max_iter')
    pars.max_iter = 500;
end
if ~isfield(pars, 'tol')
    pars.tol = 0;
end

%% Data Preprocessing
[m, n, T, C] = size(B);
N = m*n;
warning('off', 'all');

%% Parameter Initialization
lambda_1 = pars.lambda_1; lambda_2 = pars.lambda_2;
L = pars.L;

t1 = pars.t1;
t2 = 1/(8*t1*lambda_1*lambda_1);
tau = t1/(1+t1*L);

tol = pars.tol; 
max_iter = pars.max_iter;
%% Variable Initialization

X_n = A'*B;
Y_n = {zeros(m-1, n, T), zeros(m, n-1, T)};
energy_vals = zeros(max_iter+1, 1);
psnr_vals = zeros(max_iter, 1);
time_vals = zeros(max_iter, 1);
rmse_vals = zeros(max_iter, 1);
if pars.verbose
    energy_vals(1) = computeEnergy(X_n, A, B, lambda_1, lambda_2);
    psnr_vals(1) = evaluatePSNR(F_gt,X_n); 
    rmse_vals(1) = RMSE(F_gt,X_n);
    time_vals(1) = 0; 
end
tt0=tic;
%% Iteration
for i = 1:max_iter
    %
    % Iteration 
    %
    % Iterate rule (Arrow-Hurwicz and so on)
    X_bar = X_n; Y_bar = Y_n;
    Y_tilde = Y_n; 
    % Descent in primal variable
    X_b = X_bar - tau*(computeGradient(A, X_bar, B)+lambda_1*Lforward(Y_tilde)); 
    X_n_next = reshape(MatrixShrinkageOpeartor(reshape(X_b, [N, T]), tau*lambda_2), [m, n, T]);
    % Descent in dual variable
    X_tilde = 2*X_n_next-X_n;
    % Update Dual Variable
    if lambda_1 ~= 0
        Y_b = cellfun(@(X, Y) X+t2*lambda_1*Y, Y_bar, Ltrans(X_tilde), 'UniformOutput', 0);
        Y_n_next = dualNormBallProjection(Y_b);
    else
        Y_n_next = Ltrans(zeros(size(X_n)));
    end
    fprintf('Iter %d, Err: %.5f, Energy: %.5f, time: %.2f, psnr %.2f\n', i, norm(X_n_next(:) - X_n(:))/norm(X_n(:)), energy_vals(i), time_vals(i), psnr_vals(i));
    if norm(X_n_next(:) - X_n(:))/norm(X_n(:))<tol
        break
    end
    
    X_n = X_n_next; Y_n = Y_n_next;
        
    time_vals(i+1)=toc(tt0);
    
    if pars.verbose 
        % calc energy value 
        energy_vals(i+1) = computeEnergy(X_n_next, A, B, lambda_1, lambda_2);
        psnr_vals(i+1) = evaluatePSNR(F_gt,X_n_next); 
        rmse_vals(i+1) = RMSE(F_gt,X_n_next);
        if pars.debug_output
            fprintf('Iter %d, Err: %.5f, Energy: %.5f, time: %.2f, psnr %.2f\n', ...
                i, norm(X_n_next(:) - X_n(:))/norm(X_n(:)), energy_vals(i+1), time_vals(i), psnr_vals(i+1));  % Need to change X_n, since it has been changed
        end
        
    end
end

X = X_n;  

end

%% Compute Energy Function
function val = computeEnergy(X, A, B, lambda_1, lambda_2)
% TODO: Modify the energy function
    [m, n, C] = size(X);
    v = A*X-B;
    val = 1/2*(norm(v(:))^2) + ... 
        lambda_1 * TVNorm(Ltrans(X)) + ... 
        lambda_2 * NuclearNorm(reshape(X, [m*n, C]));
end

function Y_out = dualNormBallProjection(Y)
    %% Anisotropic TV: Project onto the $\ell^{\infty}$ unit ball
    [P, Q] = Y{1:2}; [m, n, T] = size(P); m = m + 1;
    N = [abs(P).^2 ;zeros(1,n,T)] + [abs(Q).^2 , zeros(m, 1, T)];
    Nr = sqrt(max(abs(N), 1));
    P = P./Nr(1:m-1, :, :); Q = Q./Nr(:, 1:n-1, :);
    Y_out = {P, Q};
end

%% Gradient of \frac{1}{2} \|AX-B\|^2
%
%   Compute Gradient of \frac{1}{2}\|A.*X-B\|^2 
%

function G = computeGradient(A, X, B)
    G = A'*(A*X-B);
end



%%  TV transformation
%
%   We reference the Lforward and Ltrans function from FISTA-TV 
%
function X=Lforward(P)
    % TODO: Modify TV calculation
    [m, n, T] = size(P{1}); m = m+1;
    X=zeros(m,n,T);
    X(1:m-1,:,:)=P{1};
    X(:,1:n-1,:)=X(:,1:n-1,:)+P{2};
    X(2:m,:,:)=X(2:m,:,:)-P{1};
    X(:,2:n,:)=X(:,2:n,:)-P{2};
end

function P=Ltrans(X)
    % TODO: Modify TV transpose calculation
    [m,n,T]=size(X);
    P{1}=X(1:m-1,:,:)-X(2:m,:,:);
    P{2}=X(:,1:n-1,:)-X(:,2:n,:);
end

%% 
%
%   Tool functions for nuclear norm and nuclear soft-thresholding
%
function val = TVNorm(Y)
    [P, Q] = Y{1:2}; [m, n, T] = size(P); m = m + 1;
    N = [abs(P).^2 ;zeros(1,n,T)] + [abs(Q).^2 , zeros(m, 1, T)];
    val = sum(sqrt(abs(N(:))));
end

function val = NuclearNorm(X)
    val = norm(svd(X, 'econ'), 1);
end

function I_tilde = MatrixShrinkageOpeartor(I, lambda)
    if lambda ~= 0
        [U, S, V] = svd(I'*I, 'econ');
        S = sqrt(S);
        U = I*V*(S^-1);
        s = diag(S);
        s_tilde = sign(s).*max(0, abs(s)-lambda); 
        I_tilde = U*diag(s_tilde)*V';
    else
        I_tilde = I;
    end
end

function X_p = RangeProjection(X, l, u)
    if((l==-Inf)&&(u==Inf))
        project=@(x)x;
    elseif (isfinite(l)&&(u==Inf))
        project=@(x)(((l<x).*x)+(l*(x<=l)));
    elseif (isfinite(u)&&(l==-Inf))
        project=@(x)(((x<u).*x)+((x>=u)*u));
    elseif ((isfinite(u)&&isfinite(l))&&(l<u))
        project=@(x)(((l<x)&(x<u)).*x)+((x>=u)*u)+(l*(x<=l));
    else
        error('lower and upper bound l,u should satisfy l<u');
    end

    X_p = project(X);
end

