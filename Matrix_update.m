function [A] = Matrix_update(B)
for i = 1:size(B,1)
    for j = 1:size(B,2)
        if B(i,j) == 0
            B(i,j) = -1;
        else
            
        end
    end
end
[A] = [B];