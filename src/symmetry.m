function H = symmetry(H_raw);
n = size(H_raw, 2);
H = H_raw;  
% take the upper trangle 
for i = 1:n
    for j = i:n
      H(i,j) = H_raw(j,i);     
    end
end
end

