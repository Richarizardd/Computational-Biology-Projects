% W=randzero(W,p)
% Randomly select a fraction of the entries of matrix W and set them to 0.
% p is the ratio of the number of zeroed entries to the total number of entries.

function W=randzero(W,p)

N=length(W(:));
nz=fix(p*N); % number of entries to be set to 0
indx=randperm(N);
indx=indx(1:nz);
W(indx)=0;