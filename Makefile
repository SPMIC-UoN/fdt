include $(FSLCONFDIR)/default.mk

PROJNAME = fdt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_CPROB} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_ZLIB}

 
DLIBS = -lmeshclass -lbint -lnewimage -lutils -lmiscmaths  -lnewmat -lfslio -lniftiio -lznz -lcprob -lprob -lm -lz


DTIFIT=dtifit
PT=probtrack
FTB=find_the_biggest
PJ=proj_thresh
MED=medianfilter
ROM=reord_OM
SAUS=sausages
DIFF_PVM=diff_pvm
RV=replacevols
MDV=make_dyadic_vectors

DTIFITOBJS=dtifit.o dtifitOptions.o
PTOBJS=probtrack.o probtrackOptions.o pt_alltracts.o pt_matrix.o pt_seeds_to_targets.o pt_simple.o pt_twomasks.o pt_matrix_mesh.o
FTBOBJS=find_the_biggest.o
PJOBJS=proj_thresh.o
MEDOBJS=medianfilter.o 
ROMOBJS=reord_OM.o
SAUSOBJS=sausages.o
DIFF_PVMOBJS=diff_pvm.o diff_pvmoptions.o
RVOBJS=replacevols.o
MDVOBJS=make_dyadic_vectors.o


SCRIPTS = eddy_correct bedpost bedpost_proc bedpost_cleanup bedpost_kill_all bedpost_kill_pid zeropad bedpost_datacheck
FSCRIPTS=correct_and_average ocmr_preproc
XFILES = dtifit probtrack find_the_biggest medianfilter diff_pvm make_dyadic_vectors proj_thresh
FXFILES = reord_OM sausages replacevols


RUNTCLS = Fdt

all: ${XFILES} ${FXFILES} 

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

${ROM}:    	${ROMOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${ROMOBJS} ${DLIBS}

${SAUS}:    	${SAUSOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SAUSOBJS} ${DLIBS}

${DIFF_PVM}:    	${DIFF_PVMOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${DIFF_PVMOBJS} ${DLIBS}

${RV}:    	${RVOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${RVOBJS} ${DLIBS}

${MDV}:    	${MDVOBJS}
		   ${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MDVOBJS} ${DLIBS}












