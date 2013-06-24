#!/bin/sh

#   Copyright (C) 2012 University of Oxford
#
#   Stamatios Sotiropoulos


devel_dir=/home/stam/fsldev/fdt  # Where the binary is 
script_dir=/home/stam/Rubix_src # Where the scripts are


if [ "x$SGE_ROOT" = "x" ] ; then
    if [ -f /usr/local/share/sge/default/common/settings.sh ] ; then
	. /usr/local/share/sge/default/common/settings.sh
    elif [ -f /usr/local/sge/default/common/settings.sh ] ; then
	. /usr/local/sge/default/common/settings.sh
    fi
fi

Usage() {
    echo ""
    echo "Usage: rubix_parallel <subject_directory> [options]"
    echo ""
    echo "ORIGINAL MODE"
    echo "expects to find bvalsLR, bvecsLR, bvalsHR, bvecsHR in subject directory"
    echo "expects to find dataLR, dataHR and nodif_brain_maskLR in subject directory"
    echo "expects to find grad_devLR and grad_devHR in subject directory, if -g is set"
    echo ""
    echo "FILTERING MODE (use -f)"
    echo "expects to find bvals, bvecs in subject directory"
    echo "expects to find data and nodif_brain_mask in subject directory"
    echo "expects to find grad_dev in subject directory, if -g is set"
    echo ""
    echo "options"
    echo "-f (runs rubix on a bedpostx directory as an anisotropic filter, default off)"
    echo "-n (number of fibres per HR voxel, default 2)"
    echo "-p (number of orientation prior modes per LR voxel, default 4)"
    echo "-w (ARD weight, more weight means less secondary fibres per voxel, default 1)"
    echo "-b (burnin period, default 5000)"
    echo "-j (number of jumps, default 1250)"
    echo "-s (sample every, default 25)"
    echo "-model (1 for monoexponential, 2 for multiexponential, default 1)"
    echo "-g (consider gradient nonlinearities, default off)"
    echo ""
    echo ""
    echo "ADDITIONALLY: you can pass on rubix options onto directly rubix_parallel"
    echo " For example:  eubix_parallel <subject directory> --fsumPrior --dPrior"
    echo " Type 'rubix --help' for a list of available options "
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
nfibres=2
nmodes=4
fudge=1
burnin=5000
njumps=1250
sampleevery=25
model=1
gflag=0
filterflag=0

shift
while [ ! -z "$1" ]
do
  case "$1" in
      -f) filterflag=1;;
      -n) nfibres=$2;shift;;
      -p) nmodes=$2;shift;;
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

opts="--nf=$nfibres --nM=$nmodes --fudge=$fudge --bi=$burnin --nj=$njumps --se=$sampleevery --model=$model"
defopts="--fsumPrior --dPrior"
opts="$opts $defopts $*"

#check that all required files exist
if [ ! -d ${subjdir} ]; then
	echo "subject directory $1 not found"
	exit 1
fi

if [ ${filterflag} -eq 1 ]; then
    echo "Running RubiX in Filter Mode"
  
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


    if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif_brain_mask` -eq 0 ]; then
	echo "${subjdir}/nodif_brain_mask not found"
	exit 1
    fi

    if [ ${gflag} -eq 1 ]; then
	if [ `${FSLDIR}/bin/imtest ${subjdir}/grad_dev` -eq 0 ]; then
	    echo "${subjdir}/grad_dev not found"
	    exit 1
	fi
    fi
else
    if [ ! -e ${subjdir}/bvecsLR ]; then
	echo "${subjdir}/bvecsLR not found"
	exit 1
    fi
    if [ ! -e ${subjdir}/bvalsLR ]; then
	echo "${subjdir}/bvalsLR not found"
	exit 1
    fi
    if [ ! -e ${subjdir}/bvecsHR ]; then
	echo "${subjdir}/bvecsHR not found"
	exit 1
    fi
    if [ ! -e ${subjdir}/bvalsHR ]; then
	echo "${subjdir}/bvalsHR not found"
	exit 1
    fi
    if [ `${FSLDIR}/bin/imtest ${subjdir}/dataLR` -eq 0 ]; then
	echo "${subjdir}/dataLR not found"
	exit 1
    fi
    if [ `${FSLDIR}/bin/imtest ${subjdir}/dataHR` -eq 0 ]; then
	echo "${subjdir}/dataHR not found"
	exit 1
    fi
    if [ `${FSLDIR}/bin/imtest ${subjdir}/nodif_brain_maskLR` -eq 0 ]; then
	echo "${subjdir}/nodif_brain_maskLR not found"
	exit 1
    fi
    if [ ${gflag} -eq 1 ]; then
	if [ `${FSLDIR}/bin/imtest ${subjdir}/grad_devLR` -eq 0 ]; then
	    echo "${subjdir}/grad_devLR not found"
	    exit 1
	fi
	if [ `${FSLDIR}/bin/imtest ${subjdir}/grad_devHR` -eq 0 ]; then
	    echo "${subjdir}/grad_devHR not found"
	    exit 1
	fi
    fi
fi

echo Making rubix directory structure

mkdir -p ${subjdir}.rubiX 
mkdir -p ${subjdir}.rubiX/diff_slices
mkdir -p ${subjdir}.rubiX/logs

nslices=0
if [ ${filterflag} -eq 1 ]; then
    nslicesHR=`$FSLDIR/bin/fslval ${subjdir}/data dim3`
    nslices=$((${nslicesHR}/2))
else
    nslices=`$FSLDIR/bin/fslval ${subjdir}/dataLR dim3`
fi

echo Queuing preprocessing stage
preprocid=`$FSLDIR/bin/fsl_sub -T 60 -R 8000 -N rubix_pre -l ${subjdir}.rubiX/logs  $script_dir/rubix_preproc.sh ${subjdir} ${filterflag} ${gflag}`
echo Queuing parallel processing stage
[ -f ${subjdir}.rubiX/commands.txt ] && rm ${subjdir}.rubiX/commands.txt
 slice=0
 while [ $slice -lt $nslices ]
 do
    slicezp=`$FSLDIR/bin/zeropad $slice 4`
    if [ ${gflag} -eq 1 ]; then
	opts="$opts --gLR=${subjdir}.rubiX/grad_devLR_slice_${slicezp} --gHR=${subjdir}.rubiX/grad_devHR_newslice_${slicezp}"
    fi    
    echo "$devel_dir/rubix --dLR=${subjdir}.rubiX/dataLR_slice_${slicezp} --mLR=${subjdir}.rubiX/nodif_brain_maskLR_slice_${slicezp} --rLR=${subjdir}.rubiX/bvecsLR --bLR=${subjdir}.rubiX/bvalsLR --dHR=${subjdir}.rubiX/dataHR_newslice_${slicezp} --rHR=${subjdir}.rubiX/bvecsHR --bHR=${subjdir}.rubiX/bvalsHR --forcedir --nf=$nfibres --nM=$nmodes --bi=$burnin --model=$model --logdir=${subjdir}.rubiX/diff_slices/data_slice_$slicezp ${opts}" >> ${subjdir}.rubiX/commands.txt
    slice=$(($slice + 1))
 done
    
 RubiXid=`${FSLDIR}/bin/fsl_sub -T 1400 -j $preprocid -l ${subjdir}.rubiX/logs -N rubix -t ${subjdir}.rubiX/commands.txt`
    
 echo Queuing post processing stage
 mergeid=`${FSLDIR}/bin/fsl_sub -T 60 -R 8000 -j $RubiXid -N rubix_post -l ${subjdir}.rubiX/logs $script_dir/rubix_postproc.sh ${subjdir} ${filterflag}`
