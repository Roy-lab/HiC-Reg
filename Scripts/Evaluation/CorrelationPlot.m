% Compute distance-stratified correlation and AUC:
% input arguments: cell, chr1, path, outpath

% cell='Gm12878'
% chr1=17
% path='.'
% outpath='.'

% bin resolution of Hi-C data:
bin=5000
% max distance in kb
distCF=1000
x=[0:bin:(1000*distCF)]/1e6;

infile=sprintf('%s/testset_error.txt',path)
% compute distance-stratified correlation, global correlation and MSE:
[ccbinRF1,CCmat1,MSEmat1]=RFdistbinEval(infile,distCF,cell,sprintf('chr%d',chr1),bin)

% compute AUC:
AUC=trapz(ccbinRF1)/(size(ccbinRF1,2)-1)
fprintf('AUC: %d\n',AUC)

% make distance-stratified correlation plot:
font=10;
pz=4;
f=figure;
L=0.5
M=0.5
plot(x,ccbinRF1,'r-o','MarkerSize',M,'LineWidth',L)  
grid on 
ylim([-0.2 1]);
box off
axis square
ylabel('Correlation','FontSize',font);
xlabel('Distance (Mbp)','FontSize',font);  
title(sprintf('%s on chr%d at %dkb',cell,chr1,bin/1000),'FontSize',font);   
set(gcf,'PaperPosition',[ 0 0 pz pz], 'PaperPositionMode','manual', 'PaperSize',[pz pz]);
saveas(gcf,sprintf('%s/HiC-Reg_%s_chr%d_%dkb_correlationplot.pdf',outpath,cell,chr1,bin/1000),'pdf');
