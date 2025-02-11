function [partition,D,overallClusters,subjectClusters] = run_kmeans(concTS,numClusters,distanceMethod,nreps,maxIter,subject_index,tsversion)

[partition,~,~,D] = kmeans(concTS,numClusters,'Distance', distanceMethod,'Replicates',nreps,'MaxIter',maxIter);

overallClusters = GET_CENTROIDS(concTS,partition,numClusters);

subjects = unique(subject_index);
nsub = length(subjects);
nparc = size(concTS,2);
subjectClusters = NaN(nsub,nparc,numClusters);

for i=1:nsub
    subjectClusters(i,:,:) = GET_CENTROIDS(concTS(subject_index==subjects(i),:),partition(subject_index==subjects(i)),numClusters);
end

if ~isempty(tsversion)
    save(fullfile([tsversion,'_k',num2str(numClusters),'_clusters.mat']),'partition','numClusters','overallClusters','subjectClusters');
end 

end