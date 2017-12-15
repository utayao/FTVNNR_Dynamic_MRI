function [ Mask ] = Gene_Mask_Matrice( F0, Mask_line )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% load perdata;
[m,n,seqlen] = size(F0);
% load Q;
line = Mask_line(1);
line1 = Mask_line(2);
N = [m,n];

% prepare sampling

OMEGA1=fftshift(MRImask(m,line1));  % set 50--0.21
k1 = length(find(OMEGA1~=0));
%k1/(m*n)
Mask1 = RandMask_InverseTransfer(OMEGA1,m,n);

% d = 5;  % other frames
% [OMEGA2] = RandMask_rect(double(m/d),double(n/d),m,n);
OMEGA2=fftshift(MRImask(m,line));  % set 50--0.21
k2 = length(find(OMEGA2~=0));
%k2/(m*n)
Mask2 = RandMask_InverseTransfer(OMEGA2,m,n);

Mask(:,:,1) = Mask1;
for i = 2:seqlen
    Mask(:,:,i) = Mask2;
end

end

