function [cellout] = cellit(func,it1,it2)
% INPUT:
% - func is function which we iterate through
%     each iteration produces a cell in our output
% - it1 is the first iterative variable and will produce
%     a cell type with one level of vars
% - it2 provides a second level of iterations; we may not need it, even
% 
% Note that to use it2, your func must have two args 
%
% OUTPUT:
% - cellout is our cell array

if nargin == 2
    cellout = cell(1,length(it1));
    Q = 1;
    for k=it1
        cellout{Q} = func(k);
        Q = Q + 1;
    end
elseif nargin == 3
    cellout = cell(length(it1),length(it2));
    Q = 1;
    for k=it1
        M = 1;
        for j=it2
            cellout{Q,M} = func(k,j);
            M = M + 1;
        end
        Q = Q + 1;
    end
else
    error('You have too few, or too many, args')
end



end

