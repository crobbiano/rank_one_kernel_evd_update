%% omega.m
function [o] = omega(x, sigma, z2, L)
%     Objective function

    o =  1 + (sigma * sum(z2(:).' ./ (L(:).' - x(:).')));