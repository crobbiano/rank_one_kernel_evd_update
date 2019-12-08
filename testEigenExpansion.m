%% Test eigen expansion stuff
clear all
addpath('expandEigen\')
%%
load('magic.mat')
X = magicgammatelescope';
datasize = 1000;
X = X(:,1:datasize);
ker = @(x,y,sigma) exp(-(sum((x-y).^2))/sigma^2);
%% Do some data things
dists = triu(pdist2(X',X'));
smean = median(dists(dists~=0));
numSampSmall = 5;
w=zeros(numSampSmall,numSampSmall);
for i=1:numSampSmall
    for j=1:numSampSmall
        w(i,j) = ker(X(:,i), X(:,j), smean);
    end
end
c = numSampSmall;
k = numSampSmall;%10;

W=w(1:c-1,1:c-1);
Wn = w;

[vW, lW] = eig(W,'vec');
[vWn, lWn] = eig(Wn,'vec');
%% Test new parts
for i=1:numSampSmall
        k(i) = ker(X(:,numSampSmall), X(:,i), smean);
end
e = k(end)/4; % k(x_new, x_new)/4;
sigma = 4/k(end);
k1 = k'; k1(end) = k1(end)/2;
k0 = k'; k0(end) = k0(end)/4;


[L, U] = expandEigensystem(lW, vW, k0(end));

[L, idxs] = sort(L);
U = U(:,idxs);
U = U(idxs,:);
k1 = k1(idxs);
k0 = k0(idxs);

[Ln,Un] = updateEigensystem(L, U, k1, sigma);
[lWnew,vWnew] = updateEigensystem(Ln, Un, k0, -sigma);
[~,ogIdxs] = sort(idxs);
vWnew = vWnew(ogIdxs,:);

%% Compare the eigenpairs calculated through EVD and rank-one update
display(['Frobenius norm of residuals between eigenvalues : ' num2str(norm(lWn - lWnew, 'fro'))])
display(['Frobenius norm of residuals between eigenvector matrices: ' num2str(norm(abs(vWn) - abs(vWnew), 'fro'))])

Wnew = vWnew*diag(lWnew)*vWnew';
display(['Frobenius norm of residuals between kernel matrices: ' num2str(norm(Wn - Wnew, 'fro'))])