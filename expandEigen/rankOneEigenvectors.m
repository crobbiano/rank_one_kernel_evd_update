%% rankOneEigenvectors.m
function [U_tilde] = rankOneEigenvectors(L, L_tilde, U, z)
%     Adjust eigenvectors
n = length(L);
LL = repmat(L,1,n);
%     LL[:,:] = np.expand_dims(L,1)
Ldiff = LL - L_tilde(:).';
LL = Ldiff.^-1;
Dz = LL.*z;
U_tilde = normc(U*Dz);

