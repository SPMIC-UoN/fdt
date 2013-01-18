#!/bin/sh

#   Copyright (C) 2012 University of Oxford
#
#   Stamatios Sotiropoulos

if [ "x$SGE_ROOT" = "x" ] ; then
    if [ -f /usr/local/share/sge/default/common/settings.sh ] ; then
	. /usr/local/share/sge/default/common/settings.sh
    elif [ -f /usr/local/sge/default/common/settings.sh ] ; then
	. /usr/local/sge/default/common/settings.sh
    fi
fi

Usage() {
    echo ""
    echo "Usage: rubix_parallel dataLR maskLR bvecsLR bvalsLR dataHR bvecsHR bvalsHR out_dir nfib nModes modelnum other_rubix_params"
    echo "       Call rubix help menu, for other_rubix_params (usually need --fsumPrior --dPrior)"
    echo ""
    exit 1
}

[ "$1" = "" ] && Usage



devel_dir=/home/stam/fsldev/fdt  # Where the binary is 
script_dir=/home/stam/Rubix_src # Where the scripts are

#option arguments
burnin=5000

dataLR=`${FSLDIR}/bin/imglob $1`; shift;  #Remove ending if exists
maskLR=`${FSLDIR}/bin/imglob $1`; shift;
bvecsLR=$1; shift;
bvalsLR=$1; shift;
dataHR=`${FSLDIR}/bin/imglob $1`; shift;
bvecsHR=$1; shift;
bvalsHR=$1; shift;
out_dir=$1; shift;
nfibres=$1; shift;
nModes=$1; shift;
modelnum=$1; shift;
opts="$*";

#check that all required files exist

if [ ! -e $bvecsLR ]; then
	echo "bvecsLR not found"
	exit 1
fi

if [ ! -e $bvalsLR ]; then
	echo "$bvalsLR not found"
	exit 1
fi


if [ ! -e $bvecsHR ]; then
	echo "bvecsHR not found"
	exit 1
fi

if [ ! -e $bvalsHR ]; then
	echo "$bvalsHR not found"
	exit 1
fi


if [ `${FSLDIR}/bin/imtest $dataLR` -eq 0 ]; then
	echo "dataLR not found"
	exit 1
fi

if [ `${FSLDIR}/bin/imtest $dataHR` -eq 0 ]; then
	echo "dataHR not found"
	exit 1
fi

if [ `${FSLDIR}/bin/imtest $maskLR` -eq 0 ]; then
	echo "brain maskLR not found"
	exit 1
fi


mkdir -p $out_dir 
mkdir -p $out_dir/diff_slices
mkdir -p $out_dir/logs


nslices=`$FSLDIR/bin/fslval $dataLR dim3`


echo Queuing preprocessing stage
preprocid=`$FSLDIR/bin/fsl_sub -T 60 -R 8000 -N rubix_pre -l $out_dir/logs  $script_dir/rubix_preproc.sh $dataLR $maskLR $dataHR $out_dir`

echo Queuing parallel processing stage
[ -f $out_dir/commands.txt ] && rm $out_dir/commands.txt
 slice=0
 while [ $slice -lt $nslices ]
 do
    slicezp=`$FSLDIR/bin/zeropad $slice 4`
    echo "$devel_dir/rubix --dLR=${out_dir}/${dataLR}_slice_${slicezp} --mLR=${out_dir}/${maskLR}_slice_${slicezp} --rLR=$bvecsLR --bLR=$bvalsLR --dHR=${out_dir}/${dataHR}_newslice_${slicezp} --rHR=$bvecsHR --bHR=$bvalsHR --forcedir --nf=$nfibres --nM=$nModes --bi=$burnin --model=$modelnum --logdir=$out_dir/diff_slices/data_slice_$slicezp ${opts}" >> $out_dir/commands.txt
    slice=$(($slice + 1))
 done
    
 RubiXid=`${FSLDIR}/bin/fsl_sub -T 1400 -j $preprocid -l $out_dir/logs -N rubix -t $out_dir/commands.txt`
    
 echo Queuing post processing stage
 mergeid=`${FSLDIR}/bin/fsl_sub -T 60 -R 8000 -j $RubiXid -N rubix_post -l $out_dir/logs $script_dir/rubix_postproc.sh $out_dir $dataLR $maskLR $dataHR`