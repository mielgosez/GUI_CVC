l = 20; %Node
neigh_i = find(AdjExt_knn(l,:) == 1); %Neighborhood of i
Ci = find(sum(VisCommunitiesExt(:,neigh_i),2) > 0); %Communities
representativity = 0;
for i = 1:length(neigh_i)
    representativity = representativity + VisCommunitiesExt(Ci(1),neigh_i(i));
end
representativity = representativity/length(neigh_i);
%% Intrinsic cohesion
ciElements = find(VisCommunitiesExt(Ci(1),:) == 1); % Elements in the first community
representativity_com = zeros(length(ciElements),1); % Representativity will be stored here

for i=1:length(ciElements)
    neigh_i = find(AdjExt_knn(ciElements(i),:) == 1); %Neighbor
    for j = 1:length(neigh_i)
        representativity_com(i) = representativity_com(i) + VisCommunitiesExt(Ci(1),neigh_i(j));
    end
    representativity_com(i) = representativity_com(i)/length(neigh_i);
end
CI = sum(representativity_com);
%% Afilliation force
ciElements = find(AdjExt_knn(l,:) == 1); % Elements in the first community
representativity_node = zeros(length(ciElements),1); % Representativity will be stored here

for i=1:length(ciElements)
    neigh_i = find(AdjExt_knn(ciElements(i),:) == 1); %Neighbor
    for j = 1:length(neigh_i)
        representativity_node(i) = representativity_node(i) + VisCommunitiesExt(Ci(1),neigh_i(j));
    end
    representativity_node(i) = representativity_node(i)/length(neigh_i);
end
FA = sum(representativity_node);