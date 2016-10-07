% Y=randsign(X,p)
% Randomly select a fraction of the entries of matrix X and flip their signs.
% p is the ratio of the number of flipped entries to the total number of entries.

function Y=randsign(X,p)

N=length(X(:));
nz=fix(p*N); % number of entries to be flipped
indx=randperm(N);
indx=indx(1:nz);
Y=X;
Y(indx)=-X(indx);