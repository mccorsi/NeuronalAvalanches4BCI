probe=activations_gr2;

for kk1=1:size(probe,2)
    for kk2=1:size(probe{kk1},2) 
        for kk3=1:size(probe{kk1}{kk2},2)
            kk4=1;
            if size(probe{kk1}{kk2}{kk3},2)>1
            sigma_temp(kk4) = sigma_avalanche(probe{kk1}{kk2}{kk3});
            kk4=kk4+1;
            else
                continue
            end
        end
        sigma_pz{kk1}(kk2)=geomean(sigma_temp);
    end
    sigma_gr(kk1)=geomean(sigma_pz{kk1});
end

%%
function [c]=sigma_avalanche(a) 
    b=sum(a,1);
    c=geomean(b(2:end)./b(1:end-1));
end
