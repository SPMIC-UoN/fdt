#!/bin/sh

subjdir=$1
filterflag=$2
gflag=$3

out_dir=${subjdir}.rubiX
maskLR=nodif_brain_maskLR

if [ ${filterflag} -eq 1 ]; then
    cp ${subjdir}/bvals $out_dir/bvalsLR
    cp ${subjdir}/bvecs $out_dir/bvecsLR
    cp ${subjdir}/bvecs $out_dir/bvecsHR
    cp ${subjdir}/bvals $out_dir/bvalsHR
    ${FSLDIR}/bin/fslslice $subjdir/data $out_dir/$dataHR
    HighRes=`fslval ${subjdir}/data pixdim1`
    LowRes=`echo "${HighRes} * 2" | bc -l`
    echo "Create Downsampled version of data at ${LowRes} mm isotropic"
    flirt -in ${subjdir}/nodif_brain_mask -ref ${subjdir}/nodif_brain_mask -applyisoxfm ${LowRes} -interp nearestneighbour -out ${out_dir}/nodif_brain_maskLR
    flirt -in ${subjdir}/data -ref ${subjdir}/data -applyisoxfm ${LowRes} -interp sinc -out ${out_dir}/dataLR
    if [ ${gflag} -eq 1 ]; then
	echo "Create Downsampled version of grad_dev at ${LowRes} mm isotropic"
	flirt -in ${subjdir}/grad_dev -ref ${subjdir}/grad_dev -applyisoxfm ${LowRes} -interp nearestneighbour -out ${out_dir}/grad_devLR
	${FSLDIR}/bin/fslslice $subjdir/grad_dev $out_dir/grad_devHR
	${FSLDIR}/bin/fslslice ${out_dir}/grad_devLR $out_dir/grad_devLR
    fi
    dataHRname=${subjdir}/data
    dataLRname=${out_dir}/dataLR
else
    cp ${subjdir}/bvalsLR $out_dir/bvalsLR
    cp ${subjdir}/bvecsLR $out_dir/bvecsLR
    cp ${subjdir}/bvecsHR $out_dir/bvecsHR
    cp ${subjdir}/bvalsHR $out_dir/bvalsHR
    ${FSLDIR}/bin/fslslice $subjdir/dataLR $out_dir/$dataLR
    ${FSLDIR}/bin/fslslice $subjdir/$maskLR $out_dir/$maskLR
    ${FSLDIR}/bin/fslslice $subjdir/dataHR $out_dir/$dataHR
    if [ ${gflag} -eq 1 ]; then
	${FSLDIR}/bin/fslslice $subjdir/grad_devLR $out_dir/grad_devLR
	${FSLDIR}/bin/fslslice $subjdir/grad_devHR $out_dir/grad_devHR
    fi
    dataLRname=${subjdir}/dataLR
    dataHRname=${subjdir}/dataHR
fi

pixdim3LR=`$FSLDIR/bin/fslval ${dataLRname} pixdim3`
pixdim3HR=`$FSLDIR/bin/fslval ${dataHRname} pixdim3`
nslicesLR=`$FSLDIR/bin/fslval ${dataLRname} dim3`
nslicesHR=`$FSLDIR/bin/fslval ${dataHRname} dim3`
zratio=`echo "scale=0; $pixdim3LR / $pixdim3HR" | bc -l`

slice=0
while [ $slice -lt $nslicesLR ]
do
    slicezp=`$FSLDIR/bin/zeropad $slice 4`
    count=$zratio
    arr_cnt=0
    while [ $count -ge 1 ]
    do 
	HRslice=`echo "( ( $slice + 1 ) * $zratio ) - $count" | bc -l`
	HRslicezp=`$FSLDIR/bin/zeropad $HRslice 4`
	filenames[arr_cnt]=${out_dir}/dataHR_slice_${HRslicezp}
	gfilenames[arr_cnt]=${out_dir}/grad_devHR_slice_${HRslicezp}
	count=$(($count - 1))
	arr_cnt=$(($arr_cnt + 1))
    done
    fslmerge -z ${out_dir}/dataHR_newslice_${slicezp} `echo ${filenames[*]}` 
    if [ ${gflag} -eq 1 ]; then
	fslmerge -z ${out_dir}/grad_devHR_newslice_${slicezp} `echo ${gfilenames[*]}` 
    fi
    slice=$(($slice + 1))
done

rm -f ${out_dir}/dataHR_slice_*
#rm -f ${out_dir}/dataLR
if [ ${gflag} -eq 1 ]; then
    rm -f ${out_dir}/grad_devHR_slice_*
#   rm -f ${out_dir}/grad_devLR
fi
echo Done
