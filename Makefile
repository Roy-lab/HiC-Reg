SRC=Framework.C FGEditor.C FGMaximizer.C GenanatomyConverter.C Move.C common/FactorManager.C common/VariableManager.C common/EvidenceManager.C common/PotentialManager.C common/FactorGraph.C common/Graph.C common/Vertex.C common/Potential.C common/RegressionTree.C common/Variable.C common/LatticeStructure.C common/SlimFactor.C common/Rule.C common/Evidence.C common/Error.C common/Distance.C
INCLPATH1=common
LIBPATH = /mnt/ws/sysbio/roygroup/shared/thirdparty/gsl_install/lib
INCLPATH2 =/mnt/ws/sysbio/roygroup/shared/thirdparty/gsl_install/include

LOCLIB=/home/dchasman/cmint_ashton_v3/learntrees/execs/learnpertarget/lib
LOCINCL=/home/dchasman/cmint_ashton_v3/learntrees/execs/learnpertarget/include

CC=g++
CFLAGS = -g -std=c++0x
LFLAG = -lgsl -lgslcblas 


#local: $(SRC)
#	$(CC) $(SRC) -I $(INCLPATH1) -I $(LOCINCL)  -L $(LOCLIB) $(LFLAG) $(CFLAGS) -o regTreeDC

ALL: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH1) -I $(INCLPATH2)  -L $(LIBPATH) $(LFLAG) $(CFLAGS) -o regForest
clean:
	rm regForest
