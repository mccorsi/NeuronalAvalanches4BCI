function [out]=func_transition_matrix(b)

%b is a binary matrix representing active region, structured as regions x
%times. Out is an asymmetric matrix, called transition matrix, containing
%the probability that, if not i is active, node j will be active in the
%next timestep.

% out_p=zeros(size(b,1),size(b,1));
% for kk1=1:size(b,1)
%     g_loop=zeros(1,size(b,1));
%     for kk2=1:size(b,2)-1
%         if b(kk1,kk2)==1
%         curr_reg=repelem(b(kk1,kk2),size(b,1));
%         next_reg=b(:,kk2+1)';
%         g=curr_reg==next_reg;
%         g_loop=g_loop+g;
%         end
%     end     
%     out_p(kk1,:)=g_loop;
% end
% out=out_p/size(b,2);

out=zeros(size(b,1),size(b,1));
num_activaz_per_reg=sum(b,2);
for kk1=1:size(b,1)
    g_loop=zeros(1,size(b,1));
    for kk2=1:size(b,2)-1
        if b(kk1,kk2)==1
        curr_reg=repelem(b(kk1,kk2),size(b,1));
        next_reg=b(:,kk2+1)';
        g=curr_reg==next_reg;
        g_loop=g_loop+g;
        end
    end    
    g_loop(g_loop~=0)=g_loop(g_loop~=0)/num_activaz_per_reg(kk1);% divides the number of times region i has been followed by ativation or region j by the total number of time region i has been active (gives a probability)
    out(kk1,:)=g_loop;
end
end