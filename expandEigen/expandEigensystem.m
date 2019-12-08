%% expandEigensystem.m
function [L, U] = expandEigensystem(L,U,e)
%     Expand eigensystem with new eigenpair orthogonal to the existing ones
%     Parameters
%     ----------
%     L : numpy.ndarray, 1d
%         Current eigenvalues
%     U : numpy.ndarray, 2d
%         Current eigenvectors
%     e : float
%         New eigenvalue

m = size(L,1);
L = [L; e];
U = [U, zeros(m,1); zeros(1,m+1)];
U(m+1,m+1) = 1;
