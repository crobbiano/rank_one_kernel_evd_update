%% updateEigensystem.m
function [L_tilde, U_tilde] = updateEigensystem(L,U,v,sigma)
%     Perform rank one update to eigensystem. Requires the eigenpairs to be
%     ordered
%
%     Parameters
%     ----------
%     L : numpy.ndarray, 1d
%         Current eigenvalues
%     U : numpy.ndarray, 2d
%         Current eigenvectors
%     v : numpy.ndarray, 2d
%         Column vector for rank one update
%     sigma : float
%         update coefficient
n = length(L);
z = U'*v;

L_tilde = rankOneEigenvalues(L, z, sigma);
U_tilde = rankOneEigenvectors(L, L_tilde, U, z);
    
