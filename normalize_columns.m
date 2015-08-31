function [A] = normalize_columns(A)
if 0
    for k=1:size(A,2)
        A(:,k) = A(:,k) / norm(A(:,k));
    end
else
    A = bsxfun(@times,A,1./sqrt(sum(abs(A).^2)));
end