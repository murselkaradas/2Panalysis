function [A, B] = deinterleave(I)
if ndims(I)<3
    A = I;
    B = [];
else
    A = I(:,:,1:2:end);
    B = I(:,:,2:2:end);
end
end