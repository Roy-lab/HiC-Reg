####################################################################################################
## Generate MULTI-CELL Features:
import os
import pandas as pd
import argparse

def ReadData(infile):
    print(infile)
    data = pd.read_table(infile, sep='\t')
    return(data)

def MergeData(wholedata):   
    Gm12878=wholedata[0]
    K562=wholedata[1]
    Huvec=wholedata[2]
    Hmec=wholedata[3]
    Nhek=wholedata[4]
    data = pd.merge(Gm12878, K562, on='Pair', suffixes=['_Gm12878', '_K562'], how='outer')  
    data1 = pd.merge(data, Huvec, on='Pair', how='outer')
    data2 = pd.merge(data1, Hmec, on='Pair', suffixes=['_Huvec', ''], how='outer')
    data3 = pd.merge(data2, Nhek, on='Pair', suffixes=['_Hmec', '_Nhek'], how='outer')   
    data4=data3.fillna(0.01)
    return(data4)


def ExtractFeat(feat,cell,outpath,filename):  
    cells=["Gm12878","K562","Huvec","Hmec","Nhek"]     
    ## colmn numbers to be removed
    ColRm=['Distance_'+s for s in cells]+['Count_'+s for s in cells]
    n = [j for j,x in enumerate(feat.columns.format()) if x in ColRm]
    #print('Deleted '+feat.columns.values[n])
    remain=list(set(range(feat.shape[1])) - set(n))
    feat1=feat.iloc[:,remain].copy()
    feat1['Distance']=feat['Distance_'+cell].values
    feat1['Distance']=feat1['Distance'].astype('int64')
    feat1['Count']=feat['Count_'+cell].values
    outname=outpath+filename
    print(outname)
    feat1.to_csv(outname, index = None,header=True,  sep='\t', mode='w',float_format='%g')


def main(args):
    inpath=args.inpath
    chr=args.chr
    ncv=args.ncv
    outpath=args.outpath
    cells=["Gm12878","K562","Huvec","Hmec","Nhek"]
    wholedata=[]
    for i in range(len(cells)):
        cell=cells[i]
        infile=inpath+'/'+cell+'_chr'+str(chr)+'_5kb_test.txt'   
        data=ReadData(infile)    
        wholedata.append(data)

    print("Merging feature data across cells")    
    datall=MergeData(wholedata)
    for i in range(len(cells)):
        cell=cells[i]    
        if not os.path.exists(outpath):
            os.makedirs(outpath)    
        n=datall.shape[0]
        validation_set_size=int(n/ncv)
        for fold in range(0,ncv):
            if fold==ncv-1:
                train_range=range(0,fold*validation_set_size)
                validation_range=range(fold*validation_set_size,n)
            else:
                train_range=range(0,fold*validation_set_size)+range((fold+1)*validation_set_size,n)
                validation_range=range(fold*validation_set_size,(fold+1)*validation_set_size)
    
            train_data=datall.iloc[train_range,:]
            test_data=datall.iloc[validation_range,:]
            file1=cell+'_MULTI-CELL_chr'+str(chr)+'_train'+str(fold)+'.txt'             
            file2=cell+'_MULTI-CELL_chr'+str(chr)+'_test'+str(fold)+'.txt'                   
            ExtractFeat(train_data,cell,outpath,file1)    
            ExtractFeat(test_data,cell,outpath,file2) 

                

if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--chr',

                        help='chromosome.',

                        type=int,

                        default=17)
    parser.add_argument('--ncv',

                        help='number of CV',

                        type=int,

                        default=1)  
    parser.add_argument('--inpath',

                        help='path of input feature data',

                        type=str,

                        default='') 
    parser.add_argument('--outpath',

                        help='path of output feature data',

                        type=str,

                        default='')     
    args = parser.parse_args()
    main(args)
    
####################################################################################################
