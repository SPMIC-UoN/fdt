include $(FSLCONFDIR)/default.mk

PROJNAME = fdt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_CPROB} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

 
DLIBS = -lmeshclass -lbint -lnewimage -lutils -lmiscmaths -lnewmat -lnewran -lfslio -lniftiio -lznz -lcprob -lprob -lm -lz
#DLIBS = -lbint -lnewimage -lutils -lmiscmaths  -lnewmat -lfslio -lniftiio -lznz -lcprob -lprob -lm -lz


DTIFIT=dtifit
CCOPS=ccops
PT=probtrack
PTX=probtrackx
FTB=find_the_biggest
PJ=proj_thresh
MED=medianfilter
ROM=reord_OM
SAUS=sausages
DIFF_PVM=diff_pvm
XFIBRES=xfibres
RV=replacevols
MDV=make_dyadic_vectors
FMO=fdt_matrix_ops
INDEXER=indexer
TEST=testfile
ORDVEC=reorder_dyadic_vectors
DPM=dpm


DTIFITOBJS=dtifit.o dtifitOptions.o
CCOPSOBJS=ccops.o ccopsOptions.o
PTOBJS=probtrack.o probtrackOptions.o pt_alltracts.o pt_matrix.o pt_seeds_to_targets.o pt_simple.o pt_twomasks.o pt_matrix_mesh.o
PTXOBJS=probtrackx.o probtrackxOptions.o streamlines.o ptx_simple.o ptx_seedmask.o ptx_twomasks.o ptx_nmasks.o
FTBOBJS=find_the_biggest.o
PJOBJS=proj_thresh.o
MEDOBJS=medianfilter.o 
ROMOBJS=reord_OM.o
SAUSOBJS=sausages.o
DIFF_PVMOBJS=diff_pvm.o diff_pvmoptions.o
XFIBOBJS=xfibres.o xfibresoptions.o
RVOBJS=replacevols.o
MDVOBJS=make_dyadic_vectors.o
FMOOBJS=fdt_matrix_ops.o
INDEXEROBJS=indexer.o
TESTOBJS=testfile.o
ORDVECOBJS=reorder_dyadic_vectors.o heap.o
DPMOBJS=dpm.o dpmOptions.o

SGEBEDPOST = bedpost bedpost_datacheck
SGEBEDPOSTX = bedpostX bedpostX_postproc.sh bedpostX_preproc.sh bedpostX_single_slice.sh

SCRIPTS = eddy_correct zeropad ${SGEBEDPOST} ${SGEBEDPOSTX}
FSCRIPTS = correct_and_average ocmr_preproc

XFILES = dtifit ccops probtrack find_the_biggest medianfilter diff_pvm make_dyadic_vectors proj_thresh
FXFILES = reord_OM sausages replacevols fdt_matrix_ops probtrackx xfibres indexer


RUNTCLS = Fdt Fdtx

all: ${XFILES} ${FXFILES} 

${PTX}:		   ${PTXOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PTXOBJS} ${DLIBS}

${PT}:		   ${PTOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PTOBJS} ${DLIBS} 

${FTB}:    	${FTBOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${FTBOBJS} ${DLIBS} 

${PJ}:    	${PJOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PJOBJS} ${DLIBS} 

${MED}:    	${MEDOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MEDOBJS} ${DLIBS} 

${DTIFIT}:    	${DTIFITOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${DTIFITOBJS} ${DLIBS}

${CCOPS}:    	${CCOPSOBJS}	
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${CCOPSOBJS} ${DLIBS}

${ROM}:    	${ROMOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${ROMOBJS} ${DLIBS}

${SAUS}:    	${SAUSOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SAUSOBJS} ${DLIBS}

${DIFF_PVM}:    	${DIFF_PVMOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${DIFF_PVMOBJS} ${DLIBS}

${XFIBRES}:    	${XFIBOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${XFIBOBJS} ${DLIBS}

${RV}:    	${RVOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RVOBJS} ${DLIBS}

${MDV}:    	${MDVOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MDVOBJS} ${DLIBS}

${FMO}:    	${FMOOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${FMOOBJS} ${DLIBS}

${INDEXER}:    	${INDEXEROBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${INDEXEROBJS} ${DLIBS}

${TEST}:    	${TESTOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${TESTOBJS} ${DLIBS}


${ORDVEC}:    	${ORDVECOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${ORDVECOBJS} ${DLIBS}

${DPM}:    	${DPMOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${DPMOBJS} ${DLIBS}












