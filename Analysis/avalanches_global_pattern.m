 
function [f,h]=avalanches_global_pattern(x,n)

% f is a matrix (containing the consecutive series of activation of one patients 
%(by consecutive I mean that all the clean bits have been attached one after the other)
% of axb, where a is the number of areas and b is the time. The argument n will specify how many
% regions to take into account. h will provide y vectors, with y the number
% of avalanches, each vector being long n, that is full with logicals where
% 0 means that that particular regions has never been active in that
% particular avalanches, and 1 means the opposite. If a region has been
% active only in 1 time bin or in more than 1 time bin, it will be 1 in the
% vector eitherway.


x=x(1:n,:);
activations_bins=find(any(x));%finds time bins with activated areas
activations_bins=[activations_bins 5 6];
gg=diff(activations_bins)~=1;
Index_start = find(gg==1);
%Index_end = strfind(gg,[1 0]);
f{1}=x(:,activations_bins(1):activations_bins(Index_start(1)));
    for ii=1:size(Index_start,2)-1 %il -1 per il 5 e 6 fittizi aggiunti a activations_bins
        f{ii+1}=x(:,activations_bins(Index_start(ii)+1):activations_bins(Index_start(ii+1)));
    end

    for kk = 1:size(f,2)
        a=sum(f{kk},2);
        h{kk}=a>0;
    end    

end