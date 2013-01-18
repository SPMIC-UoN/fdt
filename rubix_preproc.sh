#!/bin/sh

dataLR=$1
maskLR=$2
dataHR=$3
out_dir=$4

${FSLDIR}/bin/fslslice $dataLR $out_dir/$dataLR
${FSLDIR}/bin/fslslice $maskLR $out_dir/$maskLR
${FSLDIR}/bin/fslslice $dataHR $out_dir/$dataHR

pixdim3LR=`$FSLDIR/bin/fslval $dataLR pixdim3`
pixdim3HR=`$FSLDIR/bin/fslval $dataHR pixdim3`
nslicesLR=`$FSLDIR/bin/fslval $dataLR dim3`
nslicesHR=`$FSLDIR/bin/fslval $dataHR dim3`
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
	filenames[arr_cnt]=${out_dir}/${dataHR}_slice_${HRslicezp}
	count=$(($count - 1))
	arr_cnt=$(($arr_cnt + 1))
    done
    fslmerge -z ${out_dir}/${dataHR}_newslice_${slicezp} `echo ${filenames[*]}` 
    slice=$(($slice + 1))
done

rm -f ${out_dir}/${dataHR}_slice_*

echo Done
