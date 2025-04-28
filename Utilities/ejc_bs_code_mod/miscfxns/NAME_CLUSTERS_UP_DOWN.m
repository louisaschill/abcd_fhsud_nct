function [clusterNamesUp,clusterNamesDown,netangle_Up,netangle_Down] = NAME_CLUSTERS_UP_DOWN(centroids)

% Provide names for clusters based on angular distance to binary Yeo 
% Returns vector where 1 indicates a "+" state and 0 indicates a "-" state

[nparc,numClusters] = size(centroids);

if nparc == 90 %
   load('data/aal_to_yeo.csv'); 
   networklabels=aal_to_yeo;
elseif nparc == 86
    load('fs86_to_yeo.csv');
    networklabels=fs86_to_yeo;
elseif nparc == 68
    load('fs86_to_yeo.csv');
    networklabels=fs86_to_yeo(19:86);
elseif nparc ==200 || nparc ==232
   load('sch232_to_yeo.csv'); 
   networklabels=sch232_to_yeo;
elseif nparc == 100 || nparc == 116
    load('sch116_to_yeo.csv'); 
    networklabels=sch116_to_yeo;
elseif nparc == 400 || nparc == 454
    load('sch454_to_yeo.csv');
    networklabels=sch454_to_yeo;
elseif nparc == 462
    load('Lausanne_463_subnetworks.mat');
    networklabels=subnetworks;
    networklabels=subnetworks_reorder;
    networklabels(463)=[];
elseif nparc == 461
    load('Lausanne_463_subnetworks.mat');
    networklabels=subnetworks_reorder;
    networklabels([14 463])=[];
end

networklabels=networklabels(1:nparc);
networklabels=reshape(networklabels,nparc,1);

numNets = 9; %need to decide whether to include 8th network
binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == networklabels);

YeoNetNames = {'VIS', 'SOM', 'DAT', 'VAT', 'LIM', 'FPN', 'DMN','SUB', 'CER'}; %SUB = subcortical regions

% calculate cosine of angle between binary state vector and centroids
centroids_up = centroids .* (centroids > 0);
centroids_down = -1 * centroids .* (centroids < 0);     % make negative activity positive and get rid of positive activity

netangle_Up = zeros(numClusters,numNets);
netangle_Down = zeros(numClusters,numNets);

for K = 1:numClusters
    for B = 1:numNets
        netangle_Up(K,B) = dot(centroids_up(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
        netangle_Down(K,B) = dot(centroids_down(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

% get index of minimum and assign names
clusterNamesUp = cell(numClusters,1);
clusterNamesDown = cell(numClusters,1);
for K = 1:numClusters       % for up and down separately, calculate closest network
    Up_ind = find(netangle_Up(K,:) == max(netangle_Up(K,:)));
    Down_ind = find(netangle_Down(K,:) == max(netangle_Down(K,:)));    %cos(0) = 1 so need max not min (like for E.D.)
    clusterNamesUp{K} = [YeoNetNames{Up_ind},'+'];
    clusterNamesDown{K} = [YeoNetNames{Down_ind},'-'];
end