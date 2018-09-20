function [sp,flag] = checkOVLP(X,r,bac_r,overlap)

Y = X;

[idx,D] = rangesearch(X,Y,2*r); % initial determining of the overlappings
% assembling a better matrix from the one outputted by rangesearch
sp = [];
for i = 1:length(D)
    D_num = cell2mat(D(i));
    index = cell2mat(idx(i));
    for j = 2: length(D_num)
        sum_radius = bac_r(index(1)) + bac_r(index(j));
        if sum_radius - D_num(j) > overlap
            sp = [sp;index(1) index(j) ];
        end
    end
end
if ~isempty(sp)

flag = sparse(sp(:,1),sp(:,2),1,length(X(:,1)),length(X(:,1)));
else
 flag = sparse(length(X(:,1)),length(X(:,1)));   
end
end