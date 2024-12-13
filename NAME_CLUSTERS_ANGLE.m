function [clusterNames,reorderClusters,clusterNamesSort,netangle] = NAME_CLUSTERS_ANGLE(centroids)
% Provide names for clusters based on angular distance to binary Yeo
% System Vectors
% returns vector where 1 indicates a "+" state and 0 indicates a "-" state
% centroids = kClusterCentroids;

[nparc,numClusters] = size(centroids);

if nparc == 90 
    load('aal_to_yeo.csv'); 
    networklabels=aal_to_yeo;
elseif nparc == 86
    load('fs86_to_yeo.csv');
    networklabels=fs86_to_yeo;
elseif nparc == 68
     load('fs86_to_yeo.csv');
    networklabels=fs86_to_yeo(19:86);
elseif nparc == 200 || nparc == 232
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
    networklabels=subnetworks_reorder;
    networklabels(463)=[];
elseif nparc == 461 %ls463 with brainstem (last region) and region 14 removed (artefacts)
    load('data/Lausanne_463_subnetworks.mat');
    networklabels=subnetworks_reorder;
    networklabels([14 463])=[];
end

networklabels=networklabels(1:nparc);
networklabels=reshape(networklabels,nparc,1);
numNets = 9;

% make a matrix where each column corresponds to a labeled Yeo system in Lausanne parcellation
% the columns are binary vectors indicated whether a region belongs to corresponding Yeo system
binaryNetVectors = ones(nparc,numNets) .* repmat((1:numNets),[nparc 1]); 
binaryNetVectors = double(binaryNetVectors == networklabels);

% then duplicate this matrix, multiply by -1 and horizontally concatenate to
% provide separate names for when systems are low amplitude
binaryNetVectors = [binaryNetVectors, -1*binaryNetVectors];
YeoNetNames = {'VIS+', 'SOM+', 'DAT+', 'VAT+', 'LIM+', 'FPN+', 'DMN+','SUB+', 'CER+','VIS-', 'SOM-', 'DAT-', 'VAT-', 'LIM-', 'FPN-', 'DMN-','SUB-', 'CER-'};

% calculate cosine of angle between binary state vector and centroids
netangle = zeros(numClusters,numNets*2);
for K = 1:numClusters
    for B = 1:(numNets*2)
        netangle(K,B) = dot(centroids(:,K),binaryNetVectors(:,B))...
            /(norm(centroids(:,K))*norm(binaryNetVectors(:,B)));
    end
end

% get index of minimum and assign names
clusterNamesInit = cell(numClusters,1);
plusminus = true(numClusters,1);
for K = 1:numClusters
    [~, ind] = max(netangle(K,:));  
    clusterNamesInit{K} = YeoNetNames{ind};
    if ind > numNets
        plusminus(K) = false;
    end
end
clusterNames = cellstr(clusterNamesInit);

%sort by name then plus-minus
[clusterNamesSort,I] = sort(clusterNamesInit);
[~,I2] = sort(plusminus(I));
clusterNamesSort = clusterNamesSort(I2);
reorderClusters = I(I2);
