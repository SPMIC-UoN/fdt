#!/bin/sh

#   Copyright (C) 2004 University of Oxford
#
#   SHCOPYRIGHT

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FSLDIR}/lib:${CUDA}/lib64:
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH

Usage() {
    echo ""
    echo "Usage: bedpostx <subject directory> [options]"
    echo ""
    echo "expects to find bvals and bvecs in subject directory"
    echo "expects to find data and nodif_brain_mask in subject directory"
    echo "expects to find grad_dev in subject directory, if -g is set"
    echo ""
    echo "<options>:"
    echo "-NJOBS (number of jobs to queue, the data is divided in NJOBS parts, usefull for a GPU cluster, default 1)"
    echo "-n (number of fibres per voxel, default 2)"
    echo "-w (ARD weight, more weight means less secondary fibres per voxel, default 1)"
    echo "-b (burnin period, default 1000)"
    echo "-j (number of jumps, default 1250)"
    echo "-s (sample every, default 25)"
    echo "-model (Deconvolution model. 1: with sticks (default), 2: with sticks with a range of diffusivities, 3: with zeppelins)"
    echo "-g (consider gradient nonlinearities, default off)"
    echo ""
    echo ""
    echo "ALTERNATIVELY: you can pass on xfibres options onto directly bedpostx"
    echo " For example:  bedpostx <subject directory> --noard --cnonlinear"
    echo " Type 'xfibres --help' for a list of available options "
    echo " Default options will be bedpostx default (see above), and not xfibres default."
    echo ""
    echo "Note: Use EITHER old OR new syntax."
    exit 1
}

make_absolute(){
    dir=$1;
    if [ -d ${dir} ]; then
	OLDWD=`pwd`
	cd ${dir}
	dir_all=`pwd`
	cd $OLDWD
    else
	dir_all=${dir}
    fi
    echo ${dir_all}
}

[ "$1" = "" ] && Usage

subjdir=`make_absolute $1`
subjdir=`echo $subjdir | sed 's/\/$/$/g'`

echo subjectdir is $subjdir

#parse option arguments
njobs=1
nfibres=2
fudge=1
burnin=1000
njumps=1250
sampleevery=25
model=1
gflag=0

shift
while [ ! -z "$1" ]
do
  case "$1" in
      -NJOBS) njobs=$2;shift;;
      -n) nfibres=$2;shift;;
      -w) fudge=$2;shift;;
      -b) burnin=$2;shift;;
      -j) njumps=$2;shift;;
      -s) sampleevery=$2;shift;;
      -model) model=$2;shift;;
      -g) gflag=1;; 
      *) break;;
  esac
  shift
done
opts="--nf=$nfibres --fudge=$fudge --bi=$burnin --nj=$njumps --se=$sampleevery --model=$model"
defopts="--cnonlinear"
opts="$opts $defopts $*"  

#check that all required files exist

if [ ! -d $subjdir ]; then
	echo "subject directory $1 not found"
	exit 1
fi

if [ ! -e ${subjdir}/bvecs ]; then
	echo "${subjdir}/bvecs not found"
	exit 1
fi

if [ ! -e ${subjdir}/bvals ]; then
	echo "${subjdir}/bvals not found"
	exit 1
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}/data` -eq 0 ]; then
	echo "${subjdir}/data not found"
	exit 1
fi

if [ ${gflag} -eq 1 ]; then
    if [ `${FSLDIR}/bin/imtest ${subjdir}/grad_dev` -eq 0 ]; then
	echo "${subjdir}/grad_dev not found"
	exit 1
    fi
fi

if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif_brain_mask` -eq 0 ]; then
	echo "${subjdir}/nodif_brain_mask not found"
	exit 1
fi

echo Making bedpostx directory structure

mkdir -p ${subjdir}.bedpostX/
mkdir -p ${subjdir}.bedpostX/diff_parts
mkdir -p ${subjdir}.bedpostX/logs
mkdir -p ${subjdir}.bedpostX/logs/logs_gpu
mkdir -p ${subjdir}.bedpostX/logs/pid_${$}
mkdir -p ${subjdir}.bedpostX/xfms

if [ ${gflag} -eq 1 ]; then
    echo "bedpostx_multigpu_LSF "${subjdir}  $opts "-g" >> ${subjdir}.bedpostX/commands.txt
else
    echo "bedpostx_multigpu_LSF "${subjdir}  $opts >> ${subjdir}.bedpostX/commands.txt
fi

echo Copying files to bedpost directory
cp ${subjdir}/bvecs ${subjdir}/bvals ${subjdir}.bedpostX
${FSLDIR}/bin/imcp ${subjdir}/nodif_brain_mask ${subjdir}.bedpostX
if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif` = 1 ] ; then
    ${FSLDIR}/bin/fslmaths ${subjdir}/nodif -mas ${subjdir}/nodif_brain_mask ${subjdir}.bedpostX/nodif_brain
fi

part=0
echo "Submitting parts (jobs) to GPUs"

post="qsub -V -l nodes=1:ppn=1:gpus=1,walltime=8:00:00 -q dque_gpu -W depend=afterok"
while [ $part -lt $njobs ];do
	partzp=`$FSLDIR/bin/zeropad $part 4`
	
	if [ ${gflag} -eq 1 ]; then
	    	gopts="$opts --gradnonlin=${subjdir}/grad_dev"
	else
	    	gopts=$opts
	fi    
	echo "#!/bin/sh" > ${subjdir}.bedpostX/xfibres.sh
	echo "#PBS -e ${subjdir}.bedpostX/logs/error$partzp" >> ${subjdir}.bedpostX/xfibres.sh
 	echo "#PBS -o ${subjdir}.bedpostX/logs/log$partzp" >> ${subjdir}.bedpostX/xfibres.sh
	echo "${FSLDIR}/bin/xfibres_gpu --data=${subjdir}/data --mask=$subjdir.bedpostX/nodif_brain_mask -b ${subjdir}.bedpostX/bvals -r ${subjdir}.bedpostX/bvecs --forcedir --logdir=$subjdir.bedpostX/diff_parts/data_part_$partzp $gopts ${subjdir} $part $njobs" >> ${subjdir}.bedpostX/xfibres.sh
	jobid=`qsub -V -l nodes=1:ppn=1:gpus=1,walltime=8:00:00 -q dque_gpu ${subjdir}.bedpostX/xfibres.sh`
    
    	#echo $jobid
     	post=$post":${jobid}"
	
	part=$(($part + 1))
done
rm ${subjdir}.bedpostX/xfibres.sh
nvox=`${FSLDIR}/bin/fslstats $subjdir.bedpostX/nodif_brain_mask -V   | cut -d ' ' -f1 `
echo "#!/bin/sh" > ${subjdir}.bedpostX/post_proc.sh
echo "#PBS -e ${subjdir}.bedpostX/logs/error_post_proc.log" >> ${subjdir}.bedpostX/post_proc.sh
echo "#PBS -o ${subjdir}.bedpostX/logs/log_post_proc.log" >> ${subjdir}.bedpostX/post_proc.sh
echo "${FSLDIR}/bin/bedpostx_postproc_gpu.sh --data=${subjdir}/data --mask=$subjdir.bedpostX/nodif_brain_mask -b ${subjdir}.bedpostX/bvals -r ${subjdir}.bedpostX/bvecs  --forcedir --logdir=$subjdir.bedpostX/diff_parts $opts $nvox $njobs ${subjdir}" >> ${subjdir}.bedpostX/post_proc.sh
post=$post" "${subjdir}.bedpostX/post_proc.sh
#echo $post
$post
rm ${subjdir}.bedpostX/post_proc.sh
echo "All parts (jobs) submitted"

finished=0
logdir=${subjdir}.bedpostX/logs

while [ $finished -eq 0 ] ; do
    nfin=0
    part=0
    while [ $part -lt $njobs ];do
        partzp=`${FSLDIR}/bin/zeropad $part 4`
        if [ -f ${subjdir}.bedpostX/diff_parts/data_part_$partzp/mean_S0samplesJ ];then
          nfin=$(($nfin + 1))
        fi
        part=$(($part + 1))
    done
    echo $nfin "parts processed of "$njobs

    if [ -f ${subjdir}.bedpostX/xfms/eye.mat ] ; then
        finished=1
        echo "All parts processed"
    fi
    sleep 60;
done

