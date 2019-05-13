function [ccbinRF,cc,mse]=RFdistbinEval(infile,distCF,cell,prefix,bin)

d=importdata(infile);
data=d.data;
distance=data(:,end);
Ytrue=data(:,1);
Ypred=data(:,2);

% compute global correlation and MSE:
[cc,mse]=assessPerformance(Ytrue,Ypred);
fprintf('Overall CV: cc=%.3f  mse=%.3f\n',cc,mse);

% compute distance-stratified correlation:
dists=0:bin:(distCF*1000);
[paircnts,bins]=histc(distance,dists);
ccbinRF=zeros(1,size(dists,2));
for i=1:size(dists,2)
    id=find(bins==i);
    if length(id)>1
    ccb=corrcoef(Ypred(id),Ytrue(id));
    ccbinRF(i)=ccb(1,2);
    end
end
