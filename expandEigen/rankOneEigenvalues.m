%% rankOneEigenalues.m
function [L_tilde] = rankOneEigenvalues(L, z, sigma)
%     Rank one update of eigenvalues of diagonal matrix
ee=1E-12;
n = length(L);
z2 = z.^2;
z2 = z2(:,1);
factor = sigma * dot(z', z);

if sigma > 0
    bounds = [L(:).', L(end)+factor];
elseif sigma < 0
    bounds = [L(1)+factor, L(:).'];
end

L_tilde = zeros(n,1);
for i = 1:n
    a = bounds(i) + ee;
    b = bounds(i+1) - ee;
    if omega(a,sigma,z2,L) * omega(b,sigma,z2,L) > 0
        L_tilde = [];
        break;
    else
        L_tilde(i) = fzero(@(x) omega(x, sigma, z2, L), [a b]);
    end
end
end