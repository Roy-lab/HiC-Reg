function [ccbinRF,cc,mse]=RFdistbinEval(infile,distCF,cell,prefix,bin)

d=importdata(infile);
data=d.data;
distance=data(:,end);
Ytrue=data(:,1);
Ypred=data(:,2);

% compute global correlation of true and predicted counts and MSE:
cc=corrcoef(Ytrue,Ypred);
cc=cc(1,2);
if(size(Ypred,1)==1)
  mse=(sum((Ypred-Ytrue).^2))/size(Ypred,2);
else
  mse=(sum((Ypred-Ytrue).^2))/size(Ypred,1);
end

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
