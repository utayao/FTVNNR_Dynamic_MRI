function showall(im)

[m,n,T] = size(im);
im = normlize(im);
im = reshape(im,[m,n,1,T]);

figure;
montage(abs(im));