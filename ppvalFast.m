function y_out = ppvalFast(pp,xin)
% This function is a barely modified version of an answer
% to this Stack Exchange question:
% https://stackoverflow.com/questions/18642585/efficient-replacement-for-ppval
% 
% This helps remediate runtime issues because of ppval() JIT incompatibility

brk = pp.breaks.';
coef = pp.coefs;

[~, inds] = histc(xin, [-inf; brk(2:end-1); +inf]); 

x_shf = xin - brk(inds);    
zero  = ones(size(x_shf));
one   = x_shf;
two   = one .* x_shf;
three = two .* x_shf;
y_out = sum( [three two one zero] .* coef(inds,:), 2);

end
