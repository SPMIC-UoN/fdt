include $(FSLCONFDIR)/default.mk

PROJNAME = dtibayes

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_CPROB} -I${INC_PROB}  
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_CPROB} -L${LIB_PROB} 

 
DLIBS = -lmeshclass -lbint -lnewimage -lutils -lmiscmaths  -lnewmat -lavwio -lcprob -lprob -lm -lz


DTIFIT=dtifit
PT=probtrack
FTB=find_the_biggest
PJ=proj_thresh
MED=medianfilter
ROM=reord_OM
SAUS=sausages
DIFF_PVM=diff_pvm

DTIFITOBJS=dtifit.o dtifitOptions.o
PTOBJS=probtrack.o probtrackOptions.o
FTBOBJS=find_the_biggest.o
PJOBJS=proj_thresh.o
MEDOBJS=medianfilter.o 
ROMOBJS=reord_OM.o
SAUSOBJS=sausages.o
DIFF_PVMOBJS=diff_pvm.o diff_pvmoptions.o

SCRIPTS = eddy_correct
XFILES = dtifit probtrack find_the_biggest medianfilter diff_pvm
RUNTCLS = Fdt

all: ${XFILES} reord_OM sausages

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












