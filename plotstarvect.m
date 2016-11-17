function [H,padj]=plotstarvect(x,ydata,nullvalue,height,lineprops,tail)
    
    %Set default options
    defaultprops={'color','k','linestyle','none','marker','.','markersize',20};
    if nargin<5, lineprops=defaultprops; end
    if iscell(lineprops), lineprops={defaultprops{:} lineprops{:}}; end
    
    if nargin<6, tail='right'; end
    
    
    p = arrayfun(@(i) signrank(ydata(:,i),nullvalue,'tail',tail), 1:size(ydata,2));
    [~,~,padj] = fdr(p);
    p=padj<0.05;
    
    H = plot(x(p),height+0*x(p),lineprops{:});