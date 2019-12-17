function [D,P] = Matrix_update2(C,N)

for i = 1:N
    for j = 1:size(C,2)
        D(i,j) = prod(C(i+N),j:size(C,2));
        P(i,j) = prod(C(i,j:size(C,2)));
    end
end