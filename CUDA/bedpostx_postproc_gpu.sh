#!/bin/sh

#   Copyright (C) 2012 University of Oxford
#
#   SHCOPYRIGHT

subjdir=${@:${#@}}

for arg in $@;
do
	if [[ $arg =~ "--nf=" ]]; then
    		numfib=`echo $arg | cut -d '=' -f2`
	fi
done

${FSLDEVDIR}/bin/merge_parts_gpu $@

fib=1
while [ $fib -le $numfib ]
do
    ${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/merged_th${fib}samples -Tmean ${subjdir}.bedpostX/mean_th${fib}samples
    ${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/merged_ph${fib}samples -Tmean ${subjdir}.bedpostX/mean_ph${fib}samples
    ${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/merged_f${fib}samples -Tmean ${subjdir}.bedpostX/mean_f${fib}samples

    ${FSLDIR}/bin/make_dyadic_vectors ${subjdir}.bedpostX/merged_th${fib}samples ${subjdir}.bedpostX/merged_ph${fib}samples ${subjdir}/nodif_brain_mask ${subjdir}.bedpostX/dyads${fib}
    if [ $fib -ge 2 ];then
	${FSLDIR}/bin/maskdyads ${subjdir}.bedpostX/dyads${fib} ${subjdir}.bedpostX/mean_f${fib}samples
    fi

    fib=$(($fib + 1))

done

if [ `${FSLDIR}/bin/imtest ${subjdir}.bedpostX/mean_f1samples` -eq 1 ];then
    ${FSLDIR}/bin/fslmaths ${subjdir}.bedpostX/mean_f1samples -mul 0 ${subjdir}.bedpostX/mean_fsamples
    fib=1
    while [ $fib -le $numfib ]
    do
	fslmaths ${subjdir}.bedpostX/mean_fsamples -add ${subjdir}.bedpostX/mean_f${fib}samples ${subjdir}.bedpostX/mean_fsamples
	fib=$(($fib + 1))
    done	
fi



echo Removing intermediate files

if [ `imtest ${subjdir}.bedpostX/merged_th1samples` -eq 1 ];then
  if [ `imtest ${subjdir}.bedpostX/merged_ph1samples` -eq 1 ];then
    if [ `imtest ${subjdir}.bedpostX/merged_f1samples` -eq 1 ];then
      rm -rf ${subjdir}.bedpostX/diff_parts
      if [ `imtest ${subjdir}.bedpostX/grad_dev_slice_0000` -eq 1 ];then
	  rm -f ${subjdir}.bedpostX/grad_dev_slice_*
      fi	  
    fi
  fi
fi

echo Creating identity xfm

xfmdir=${subjdir}.bedpostX/xfms
echo 1 0 0 0 > ${xfmdir}/eye.mat
echo 0 1 0 0 >> ${xfmdir}/eye.mat
echo 0 0 1 0 >> ${xfmdir}/eye.mat
echo 0 0 0 1 >> ${xfmdir}/eye.mat

echo Done
