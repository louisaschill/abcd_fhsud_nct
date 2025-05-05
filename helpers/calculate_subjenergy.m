%% Compute subject-specific energy matrices after choosing T with T_sweep_sps.m
function [E_full,E_region, globalTE, regionalTE] = calculate_subjenergy(SC, sub_centroids,subj_scanInd, c, T, numClusters)
% energy matrices: 
% define x0 and xf, initial and final states as cluster centroids for each state transition

if size(SC,1) == size(SC,2) % AVERAGE SC 
    nsubjs = size(sub_centroids,1);
    nparc = size(sub_centroids, 2);

    Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
    Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
    onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
    offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

    E_full = NaN(nsubjs,numClusters^2);
    E_region = NaN(nsubjs,numClusters^2,nparc);

    for i = 1:nsubjs
        SC(1:nparc+1:end) = 0; % DIAGONOL ZEROS
        Anorm= NORMALIZE(SC,c);
        sub_centroid = squeeze(sub_centroids(i,:,:));
        x0 = squeeze(sub_centroid(:,Xo_ind));
        xf = squeeze(sub_centroid(:,Xf_ind)); % now each column of x0 and xf represent state transitions
        WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
        [E_full(i,:),E_region(i,:,:)] = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each
    end

else % INDIVIDUAL SC 
    nsubjs = size(sub_centroids,1);
    nparc = size(sub_centroids, 2);

    if size(SC{1},1) ~= nparc
        nparc = 68; 
        sub_centroids = sub_centroids(:,19:86,:);
    end

    Xf_ind = repmat(1:numClusters,[1 numClusters]); % final state order
    Xo_ind = repelem(1:numClusters,numClusters); % paired with different initial states, use reshape(x,[numClusters numClusters])' to get matrix
    onDiag = (1:numClusters) + (numClusters*(0:(numClusters-1)));
    offDiag = 1:(numClusters^2); offDiag(onDiag) = []; % isolate off diagonal from linearized transition probabilities

    E_full = NaN(nsubjs,numClusters^2);
    E_region = NaN(nsubjs,numClusters^2,nparc);

    for i = 1:nsubjs
        sc = SC{i};
        sc(1:nparc+1:end) = 0; % DIAGONOL ZEROS
        Anorm= NORMALIZE(sc,c);
        sub_centroid = squeeze(sub_centroids(i,:,:));
        x0 = squeeze(sub_centroid(:,Xo_ind));
        xf = squeeze(sub_centroid(:,Xf_ind)); % now each column of x0 and xf represent state transitions
        WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T
        [E_full(i,:),E_region(i,:,:)] = MIN_CONTROL_ENERGY(Anorm, WcI, x0, xf, T,false); % compute minimum control energy for each
    end 
end 
globalTE = mean(E_full,2); 
regionalTE = squeeze(mean(E_region, 2));
