function out = normlize(x)
x = abs(x);

[m,n,T] = size(x);

for t=1:T
tp = x(:,:,t);
tp = tp(:);
% u = mean(tp);
% s = std(tp);
l = min(tp);
u = max(tp);

y = (tp-l)/(u-l);
out(:,:,t) =reshape(y,[m,n]);


end