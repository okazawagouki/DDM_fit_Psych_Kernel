
function newd = local_sum ( oldd , sumlen , Dim)

% function  B = local_sum(A, L, dim)
%
% B=local_sum(A,L) breaks vector A into chuncks of length L and claulates the sum elements
% in each chunck. The length of B will be 1/L of the length of A. L must be an integer.
%
% if A is a matrix, the operation is done along the dimension specified by
% dim (default 1).
%

if isvector(oldd)
    newd_len = floor(length(oldd)/sumlen);
    oldd_len = newd_len * sumlen;

    newd = sum(reshape(oldd(1:oldd_len),[sumlen newd_len]))';
else
    % if oldd is matrix, then local sum along Dim
    if ~exist('Dim', 'var')
        Dim = 1;
    end
    if Dim > ndims(oldd)
        error('dim is greater than the number of dimensions of input matrix A');
    end
    didx = 1:ndims(oldd);
    didx([1, Dim]) = didx([Dim, 1]); 
    oldd = permute(oldd, didx); % swap target dimension to 1
    
    newd_len = floor(size(oldd,1)/sumlen);
    if newd_len == 0
        error('The size of the input matrix A is smaller than L');
    end
    oldd_len = newd_len * sumlen;
    full_size = size(oldd);
    newd = sum(reshape(oldd(1:oldd_len, :), [sumlen, newd_len, full_size(2:end)]), 1);
    newd = reshape(newd, [newd_len, full_size(2:end)]);
    newd = permute(newd, didx); % swap back
end


