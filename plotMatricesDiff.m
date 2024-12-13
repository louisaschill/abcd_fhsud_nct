function [f,grpAvg1, grpAvg2] = plotMatricesDiff(data, groupingVariable,numRows,...
    rowNames,saveName,nameOrder,titletext,p,bhfdr,colorLim,clusterOrder,diffOnly,colorbarTitle)
% most be comparing two categories only

if isempty(clusterOrder)
    clusterOrder = 1:numRows;
end
if ~isempty(nameOrder)
    name1 = categorical(nameOrder(1));
    name2 = categorical(nameOrder(2));
else
    uniqueCategories = unique(groupingVariable,'stable');
    name1 = uniqueCategories(1);
    name2 = uniqueCategories(2);
end

data1 = data(groupingVariable==name1,:);
data2 = data(groupingVariable==name2,:);
[~,pavg,~,t]=ttest2(data1,data2);

if isempty(p)
    fdravg = mafdr(pavg,'BHFDR',bhfdr);
    fdravg = reshape(fdravg,[numRows numRows])';
    pavg = reshape(pavg,[numRows numRows])';
    fdravg = fdravg(clusterOrder, clusterOrder);
    pavg = pavg(clusterOrder, clusterOrder);

else
    fdravg = mafdr(p,'BHFDR',bhfdr);
    fdravg = reshape(fdravg,[numRows numRows])';
    pavg = reshape(p,[numRows numRows])';
    fdravg = fdravg(clusterOrder, clusterOrder);
    pavg = pavg(clusterOrder, clusterOrder);
end

grpAvg1 = reshape(nanmean(data1),[numRows numRows]);
grpAvg2 = reshape(nanmean(data2),[numRows numRows]);
grpAvg1 = grpAvg1(clusterOrder, clusterOrder);
grpAvg2 = grpAvg2(clusterOrder, clusterOrder);

grpDiff = reshape(squeeze(t.tstat),[numRows numRows])';% .* -log(fdravg); (add sign() around tstat if want
grpDiff = grpDiff(clusterOrder, clusterOrder);

maxVal = max(max([grpAvg2,grpAvg1])); % sync color scales
minVal = min(min([grpAvg2,grpAvg1]));

if ~diffOnly
    f=figure; set(f,'Position',[0,100,1000,500]);
    fontsize = 25;

    subplot(1,3,1);
    imagesc(grpAvg2);
    xticks(1:numRows); yticks(1:numRows); colormap('plasma');
    xticklabels(rowNames(clusterOrder)); xtickangle(90); yticklabels(rowNames(clusterOrder)); axis square;
    COLOR_TICK_LABELS(true,true,numRows);
    ylabel('initial state'); xlabel('final state');
    if isempty(titletext)
        title(name1,'fontsize',fontsize);
    end 
    set(gca,'FontSize',fontsize);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','calibri');
    clim([minVal maxVal]); colorbar('FontName','calibri');

    subplot(1,3,2);
    imagesc(grpAvg1);
    xticks(1:numRows); yticks(1:numRows); colormap('plasma');
    xticklabels(rowNames(clusterOrder)); xtickangle(90); yticklabels(rowNames(clusterOrder)); axis square;
    COLOR_TICK_LABELS(true,true,numRows);
    ylabel('initial state'); xlabel('final state');
    if isempty(titletext)
        title(name2,'fontsize',fontsize);
    end 
    set(gca,'FontSize',fontsize);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','calibri');
    clim([minVal maxVal]); colorbar('FontName','calibri');

    subplot(1,3,3);
    imagesc(grpDiff); colormap('plasma');
    xticks(1:numRows); xticklabels(rowNames(clusterOrder)); xtickangle(90);
    yticks(1:numRows); yticklabels(rowNames(clusterOrder)); axis square
    ylabel('initial state'); xlabel('final state');
    sig_thresh = 0.05;
    [y,x] = find(pavg < sig_thresh);
    text(x-.25,y+.2,'*','Color',[0.8,0.8,0.8],'Fontsize', 60);
    [y,x] = find(fdravg < sig_thresh);
    text(x-.25,y+.2,'**','Color',[0.8,0.8,0.8],'Fontsize', 60);
    h = colorbar('FontName','calibri'); ylabel(h,colorbarTitle);
    if isempty(colorLim)
        u_caxis_bound = max(max(grpDiff));
        l_caxis_bound = min(min(grpDiff));
        clim([l_caxis_bound u_caxis_bound]);
        h.Ticks = [l_caxis_bound (u_caxis_bound+l_caxis_bound)/2 u_caxis_bound];
        h.TickLabels = [round(l_caxis_bound,2,'significant') round((l_caxis_bound+u_caxis_bound)/2,2,'significant') round(u_caxis_bound,1,'significant')];

    else
        clim(colorLim);
        h.Ticks = [colorLim(1) (colorLim(2) +colorLim(1))/2 colorLim(2)];
        h.TickLabels = [round(colorLim(1),2,'significant') round((colorLim(1)+colorLim(2))/2,2,'significant') round(colorLim(2),1,'significant')];
    end
    COLOR_TICK_LABELS(true,true,numRows);
    
    if isempty(titletext)
        title([char(name1) ' - ' char(name2)],'fontsize',fontsize);
        set(gca,'FontSize',fontsize);
    else
        sgtitle(titletext,'FontSize',fontsize);
    end
    set(gca,'FontSize',fontsize);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','calibri');

else
    f=figure; set(f, 'Position', [100,100,700,600]);
    fontsize = 25; 

    imagesc(grpDiff); 
    xticks(1:numRows); xticklabels(rowNames(clusterOrder)); xtickangle(90);
    yticks(1:numRows); yticklabels(rowNames(clusterOrder)); axis square
    ylabel('initial state'); xlabel('final state');
    sig_thresh = 0.05;
    [y,x] = find(pavg < sig_thresh);
    text(x-.25,y+.2,'*','Color',[1 1 1],'Fontsize', 60);
    [y,x] = find(fdravg < sig_thresh);
    text(x-.25,y+.2,'**','Color',[1 1 1],'Fontsize', 60);
    h = colorbar('FontName','calibri'); ylabel(h,colorbarTitle);
    if isempty(colorLim)
        u_caxis_bound = max(max(grpDiff));
        l_caxis_bound = min(min(grpDiff));
        caxis([l_caxis_bound u_caxis_bound]);
        h.Ticks = [l_caxis_bound (u_caxis_bound+l_caxis_bound)/2 u_caxis_bound];
        h.TickLabels = [round(l_caxis_bound,3,'significant') round((l_caxis_bound+u_caxis_bound)/2,3,'significant') round(u_caxis_bound,3,'significant')];
    else
        l_caxis_bound = colorLim(1); u_caxis_bound = colorLim(2);
        clim(colorLim);
        h.Ticks = [colorLim(1) (colorLim(2) +colorLim(1))/2 colorLim(2)];
        h.TickLabels = [round(colorLim(1),3,'significant') round((colorLim(1)+colorLim(2))/2,3,'significant') round(colorLim(2),3,'significant')];
    end
    
    %colormap('parula');
    colormap(flipud(cbrewer2('RdYlBu')));
   COLOR_TICK_LABELS(true,true,numRows);
    if isempty(titletext)
        title([char(name1) ' - ' char(name2)],'fontsize',fontsize);
    else 
        title(titletext,'fontsize',fontsize);
    end 
    
    set(gca,'FontSize',fontsize);
    set(gca,'TickLength',[0 0]);
    set(gca,'Fontname','calibri');
end

set(gca,'Fontname','calibri');
set(gca,'FontSize',fontsize);

if ~isempty(saveName)
    saveas(f,saveName);
end
