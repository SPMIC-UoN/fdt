include $(FSLCONFDIR)/default.mk

PROJNAME = dtibayes

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_CPROB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_CPROB} 

 
DLIBS = -lnewimage -lutils -lmiscmaths  -lnewmat -lavwio -lcprob -lm -lz


DTIFIT=dtifit
PT=probtrack
FTB=find_the_biggest
PJ=proj_thresh
MED=medianfilter
ROM=reord_OM
SAUS=sausages

DTIFITOBJS=dtifit.o dtifitOptions.o
PTOBJS=probtrack.o probtrackOptions.o
FTBOBJS=find_the_biggest.o
PJOBJS=proj_thresh.o
MEDOBJS=medianfilter.o 
ROMOBJS=reord_OM.o
SAUSOBJS=sausages.o

SCRIPTS = eddy_correct
XFILES = dtifit probtrack find_the_biggest medianfilter
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






