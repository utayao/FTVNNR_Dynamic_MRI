function sol=RMSE(M,F)
 
% M = abs(M(:));
% F = abs(F(:));
% 
% 
% M=double(M);
% F=double(F);
% [n,m,d]=size(F);
D=(M-F); 
% sol=sqrt(sum(sum(sum(D)))/(n*m*d));

sol=sqrt(sum(abs(D(:)).^2)/sum(abs(F(:)).^2));
% err = M-F;
% err=double(err);
% [n,m]=size(err);
% D = abs(err).^2;
% sol=sqrt(sum(sum(sum(D)))/(n*m));
