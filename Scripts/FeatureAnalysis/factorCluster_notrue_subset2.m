function newd_reordered=factorCluster(cellline,featurefile,outfile,outdir)
d=importdata(featurefile);
data=d.data;
name0=strrep(d.textdata(1,2:end),'-','');
name=strrep(name0,'_','-');
nparts=5;
%[newd,u,v,err]=nmf_smoothing(data,nparts);
restarts=20;%20;
%restarts=10;
ubest=[];
vbest=[];
besterr=1000000000;
for rst=1:restarts
[u,v]=nnmf(data,nparts);
newd=u*v;
err=sum((newd-data).^2);
if(err<besterr)
ubest=u;
vbest=v;
besterr=err;
fprintf('updating with solution at %d\n',rst);
end
fprintf('restart iter=%d\n',rst);
end
u=ubest;
v=vbest;
newd=u*v;

minval=min(data(:));
zs=find(data<=minval);
%newd(zs)=minval;
u=u./repmat(sum(u,2),1,size(u,2));
v=v./repmat(sum(v,1),size(v,1),1);
%figure(1);
figure;
%now first assign a gene/cell to  a dimension and then reorder
[ig,rid]=max(u,[],2);
[ig,cid]=max(v);
rowstopick=[];
colstopick=[];
for i=1:size(u,1)
[ig,order]=sort(u(i,:),'descend');
margu(i)=u(i,order(1))-u(i,order(2));
%if(u(i,order(1))-u(i,order(2))>0.01)
if(u(i,order(1))-u(i,order(2))>0.1)
rowstopick=[rowstopick i];
end
end


for i=1:size(v,2)
[ig,order]=sort(v(:,i),'descend');
margv(i)=v(order(1),i)-v(order(2),i);
%if(v(order(1),i)-v(order(2),i)>0.02)
if(v(order(1),i)-v(order(2),i)>0.1)
colstopick=[colstopick i];
end
end

fprintf('Picked %d of %d rows and %d of %d cols\n',size(rowstopick,2),size(u,1),size(colstopick,2),size(v,2));

[ig,rid]=max(u(rowstopick,:),[],2);
[ig,cid]=max(v(:,colstopick));


[ig,rorder]=sort(rid);
[ig,corder]=sort(cid');
rparts=[];
cparts=[];
rparts=getparts(rid,rorder)
cparts=getparts(cid',corder)

h=subplot(2,2,2);
set(h,'Position',[0.2 0.8 0.75 0.1]);
%imagesc(v(:,colstopick(corder))',[0 max(v(:))]);
imagesc(v(:,colstopick(corder)),[0 max(v(:))]);
for i=1:nparts-1
%line([0, nparts+0.5],[cparts(i) cparts(i)],'color',[1 1 1]);
line([cparts(i) cparts(i)],[0 nparts+0.5],'color',[1 1 1]);
end
%imagesc(u(rorder,:),[0 1]);
h=subplot(2,2,3);
set(h,'position',[0.05 0.05 0.1 0.7])
%imagesc(v(:,corder)',[0 1]);
imagesc(u(rowstopick(rorder),:),[0 max(u(:))]);
%title('max dim');
for i=1:nparts-1
line([0 nparts+0.5],[rparts(i) rparts(i)],'color',[1 1 1]);
end

h=subplot(2,2,4);
set(h,'position',[0.2 0.05 0.75 0.7])
%imagesc(newd(rorder,corder),[0.8 1.2]);
newd_reordered=newd(rowstopick(rorder),colstopick(corder));
%imagesc(newd_reordered',[0 5]);
imagesc(newd_reordered,[0 5]);
set(gca,'xtick',1:length(colstopick),'xticklabel',name(colstopick(corder)),'FontSize',3,'XTickLabelRotation',90,'TickLength',[0 0]);

%imagesc(newd(rorder,corder),[.8 1.2]);
for i=1:nparts-1
line([0 size(newd,2)],[rparts(i) rparts(i)],'color',[1 1 1],'linewidth',1);
%line([rparts(i) rparts(i)], [0 size(newd_reordered,2)],'color',[1 1 1],'linewidth',1);
end
hold on;
for i=1:nparts-1
line([cparts(i) cparts(i)],[0 size(newd,1)],'color',[1 1 1],'linewidth',1);
%line([0 size(newd_reordered,1)],[cparts(i) cparts(i)],'color',[0 0 0],'linewidth',1);
end

saveas(gcf,sprintf('%s/nmf_out_%s.pdf',outdir,cellline),'pdf');

genenames=d.textdata(2:end,1);
featnames=d.textdata(1,2:end);
fid=fopen(sprintf('%s/%s_nmf_maxdim_cellclust.txt',outdir,cellline),'w');
for i=1:size(corder,1)
imp=sum(newd(rowstopick,colstopick(corder(i))));
fprintf(fid,'%s\t%d\t%f\n',featnames{colstopick(corder(i))},cid(corder(i)),imp);
end
fclose(fid);
fid=fopen(sprintf('%s/%s_nmf_maxdim_geneclust.txt',outdir,cellline),'w');
for i=1:size(rorder,1)
fprintf(fid,'%s\t%d\n',genenames{rowstopick(rorder(i))},rid(rorder(i)));
end
fclose(fid);
%overlayNMF_Tsne(newd,rid);
%return

%figure(2);

allrowmeans=[];
for c=1:nparts
mems=find(rid==c);
if(length(mems)==1)
    cmean=newd(rowstopick(mems),:);
else
    cmean=mean(newd(rowstopick(mems),:));
end
    allrowmeans=[allrowmeans; cmean];
end

%fid=fopen('nmfclustered_maxid_out_top5_noltgt_K562.txt','w');
fid=fopen(outfile,'w');
%fprintf(fid,'CID');
%for c=1:size(allrowmeans,1)
%fprintf(fid,' CID%d',c);
%end
%fprintf(fid,'\n');
currcid=cid(corder(1));
for c=1:size(corder,1)
for k=1:size(allrowmeans,1)
%fprintf(fid,'%s',featnames{colstopick(corder(c))});
%fprintf(fid,'%s',featnames{order(c)});
%overlayNMF_Tsne(newd,rid);
%fprintf(fid,' %.2f',allrowmeans(k,colstopick(corder(c))));
%fprintf(fid,' %.2f',newrowmeans(order(c),k));
fprintf(fid,'%s||8 CID%d||8 %.2f\n',featnames{colstopick(corder(c))},k,allrowmeans(k,colstopick(corder(c))));
end
newcid=cid(corder(c));
if(newcid~=currcid)
fprintf(fid,'|Spacer||ClusterSpacer%d\n',newcid);
currcid=newcid;
end
%fprintf(fid,'\n');
end
fclose(fid);
