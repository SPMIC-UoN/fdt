include $(FSLCONFDIR)/default.mk

PROJNAME = fdt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_CPROB} -I${INC_PROB} -I${INC_BOOST} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

DLIBS = -lwarpfns -lbasisfield -lmeshclass -lnewimage -lutils -lmiscmaths -lnewmat -lnewran -lfslio -lniftiio -lznz -lcprob -lprob -lm -lz

DTIFIT=dtifit
CCOPS=ccops
PTX=probtrackx
MED=medianfilter
ROM=reord_OM
SAUS=sausages
XFIBRES=xfibres
XFIBRES2=xfibres_2
RV=replacevols
MDV=make_dyadic_vectors
FMO=fdt_matrix_ops
INDEXER=indexer
TEST=testfile
VECREG=vecreg
KURTOSIS=kurtosis
SWAPDYADS=swap_dyadic_vectors
PVMFIT=pvmfit
DTIGEN=dtigen
RARNG=rearrange
XPRED=xfibres_pred
RUBIX=rubix


DTIFITOBJS=dtifit.o dtifitOptions.o diffmodels.o Bingham_Watson_approx.o
CCOPSOBJS=ccops.o ccopsOptions.o
PTXOBJS=probtrackx.o probtrackxOptions.o streamlines.o ptx_simple.o ptx_seedmask.o ptx_twomasks.o ptx_nmasks.o ptx_meshmask.o
MEDOBJS=medianfilter.o 
ROMOBJS=reord_OM.o
SAUSOBJS=sausages.o
XFIBOBJS=xfibres.o xfibresoptions.o diffmodels.o Bingham_Watson_approx.o 
XFIBOBJS2=xfibres_2.o xfibresoptions.o
RVOBJS=replacevols.o
MDVOBJS=make_dyadic_vectors.o
FMOOBJS=fdt_matrix_ops.o
INDEXEROBJS=indexer.o
TESTOBJS=testfile.o 
VECREGOBJS=vecreg.o
KURTOSISOBJS=kurtosis.o dtifitOptions.o
SWAPDYADSOBJS=swap_dyadic_vectors.o
PVMFITOBJS=pvmfit.o pvmfitOptions.o diffmodels.o Bingham_Watson_approx.o
DTIGENOBJS=dtigen.o
RARNGOBJS=rearrange.o
XPREDOBJS=xfibres_pred.o
RUBIXOBJS=rubix.o diffmodels.o rubixvox.o rubixoptions.o Bingham_Watson_approx.o

SGEBEDPOST = bedpost 
SGEBEDPOSTX = bedpostx bedpostx_postproc.sh bedpostx_preproc.sh bedpostx_single_slice.sh bedpostx_datacheck

SCRIPTS = eddy_correct zeropad maskdyads probtrack fdt_rotate_bvecs ${SGEBEDPOST} ${SGEBEDPOSTX} 
FSCRIPTS = correct_and_average ocmr_preproc

XFILES = dtifit ccops medianfilter make_dyadic_vectors vecreg xfibres probtrackx pvmfit dtigen

FXFILES = reord_OM sausages replacevols fdt_matrix_ops indexer rearrange xfibres_pred


RUNTCLS = Fdt

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

${XFIBRES}:    	${XFIBOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${XFIBOBJS} ${DLIBS}

${XFIBRES2}:    	${XFIBOBJS2}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${XFIBOBJS2} ${DLIBS}

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

${VECREG}:    	${VECREGOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${VECREGOBJS} ${DLIBS}


${KURTOSIS}:   ${KURTOSISOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${KURTOSISOBJS} ${DLIBS}

${SWAPDYADS}: ${SWAPDYADSOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SWAPDYADSOBJS} ${DLIBS}

${PVMFIT}:    	${PVMFITOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PVMFITOBJS} ${DLIBS}

${DTIGEN}:    	${DTIGENOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${DTIGENOBJS} ${DLIBS}

${RARNG}: 	${RARNGOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RARNGOBJS} ${DLIBS}

${XPRED}: 	${XPREDOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${XPREDOBJS} ${DLIBS}

${RUBIX}: 	${RUBIXOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RUBIXOBJS} ${DLIBS}