%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an example code of FTVNNR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you find this code is useful for your research, please cite the following papers:
%
% [1] Jiawen Yao, Zheng Xu, Xiaolei Huang, Junzhou Huang  
% " An Efficient Algorithm for Dynamic MRI Using Low-Rank and Total Variation Regularizations",
% In Medical Image Analysis, 44, 14-27, 2018.
% [2] Jiawen Yao, Zheng Xu, Xiaolei Huang, Junzhou Huang 
% "Accelerated dynamic MRI reconstruction with total variation and nuclear norm regularization", 
% In MICCAI 2015.
% 
% The code also include partial useful functions from SSMRI MATLAB Toolbox
% http://web.engr.illinois.edu/~cchen156/SSMRI.html
%
% By Jiawen Yao@UT Arlington
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clear all;   
addpath(genpath('.'))

%% Load Data
% The data is from SSMRI Toolbox 
% 
load perdata.mat
F_gt = perdata;
[m,n,T,c] = size(F_gt); 
N = [m,n]; 

%% Generate Mask
pars.image_size = N;
pars.d = 4;
for i = 1:T 
    mask(:,:,i) = load_mask( 'cartesian', pars );    
end

A = multi_p2DFT(mask);
B = A*F_gt;

showall(mask);

%% FTVNNR
tic  
pars.lambda_1 = 0.001; 
pars.lambda_2 = 2; 
xhat_TVLR = Solve_TVLR(A, B, pars, F_gt);
Time_TVLR = toc;            

% Results
% showall(F_gt);
% showall(xhat_TVLR);

% Show results
for i = 1:size(xhat_TVLR,3)
     figure(1); clf;
     subplot(1,2,1);
     imshow(normlize(F_gt(:,:,i))), axis off, colormap gray; axis off;
     title(['Original image:',num2str(i)],'fontsize',12);
     
     subplot(1,2,2);
     imshow(normlize(xhat_TVLR(:,:,i))), axis off,colormap gray; axis off;
     title('FTVNNR','fontsize',12);   
     pause(0.04)
end
