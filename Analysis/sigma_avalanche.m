
function [c]=sigma_avalanche(a) 
    b=sum(a,1);
    c=geomean(b(2:end)./b(1:end-1));
end