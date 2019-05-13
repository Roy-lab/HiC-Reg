function [cc,mse]=assessPerformance(ypred, ytrue)
% compute the correlation and mse of true and predicted counts
cc=corrcoef(ypred,ytrue);
cc=cc(1,2);
if(size(ypred,1)==1)
  mse=(sum((ypred-ytrue).^2))/size(ypred,2);
else
  mse=(sum((ypred-ytrue).^2))/size(ypred,1);
end
