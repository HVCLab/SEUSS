function [plv] = ecog_plv(angles, dim)
if nargin<2 
    dim = 1; 
end
plv = abs(nansum(exp(1i*angles),dim))/size(angles,dim);
