#!/bin/sh

subjdir=$1
filterflag=$2

out_dir=${subjdir}.rubiX
maskLR=nodif_brain_maskLR


numfib=`${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_0000/f*samples* | wc -w | awk '{print $1}'`

fib=1
while [ $fib -le $numfib ]
do
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/merged_th${fib}samples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/th${fib}samples*`
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/merged_ph${fib}samples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/ph${fib}samples*`
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/merged_f${fib}samples  `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/f${fib}samples*`
    ${FSLDIR}/bin/fslmaths ${out_dir}/merged_f${fib}samples -Tmean ${out_dir}/mean_f${fib}samples
    ${FSLDIR}/bin/make_dyadic_vectors ${out_dir}/merged_th${fib}samples ${out_dir}/merged_ph${fib}samples ${out_dir}/mean_f${fib}samples ${out_dir}/dyads${fib} 95
    if [ $fib -ge 2 ];then
	${FSLDIR}/bin/maskdyads ${out_dir}/dyads${fib} ${out_dir}/mean_f${fib}samples
    fi
    fib=$(($fib + 1))
done

if [ `${FSLDIR}/bin/imtest ${out_dir}/mean_f1samples` -eq 1 ];then
    ${FSLDIR}/bin/fslmaths ${out_dir}/mean_f1samples -mul 0 ${out_dir}/mean_fsumsamples
    fib=1
    while [ $fib -le $numfib ]
    do
	${FSLDIR}/bin/fslmaths ${out_dir}/mean_fsumsamples -add ${out_dir}/mean_f${fib}samples ${out_dir}/mean_fsumsamples
	fib=$(($fib + 1))
    done	
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/mean_dsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/mean_dsamples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/mean_dsamples*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/mean_S0samples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/mean_S0samples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/mean_S0samples*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/mean_tausamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/mean_tausamples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/mean_tausamples*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/mean_tau_LRsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/mean_tau_LRsamples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/mean_tau_LRsamples*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/mean_S0_LRsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/mean_S0_LRsamples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/mean_S0_LRsamples*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/En_Lik_LRsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/En_Lik_LRsamples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/En_Lik_LRsamples*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/En_Prior_LRsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/En_Prior_LRsamples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/En_Prior_LRsamples*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/mean_dstd_samples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/mean_dstd_samples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/mean_dstd_samples*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/HRbrain_mask` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/HRbrain_mask `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/HRbrain_mask*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/mean_fsumPriorMode_LRsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/mean_fsumPriorMode_LRsamples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/mean_fsumPriorMode_LRsamples*`
fi

if [ `${FSLDIR}/bin/imtest ${out_dir}/diff_slices/data_slice_0000/mean_dPriorMode_LRsamples` -eq 1 ];then
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/mean_dPriorMode_LRsamples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/mean_dPriorMode_LRsamples*`
fi

numModes=`${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_0000/Pk*samples* | wc -w | awk '{print $1}'`

mode=1
while [ $mode -le $numModes ]
do
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/merged_Pth${mode}samples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/Pth${mode}samples*`
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/merged_Pph${mode}samples `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/Pph${mode}samples*`
    ${FSLDIR}/bin/fslmerge -z ${out_dir}/merged_Pk${mode}samples  `${FSLDIR}/bin/imglob ${out_dir}/diff_slices/data_slice_*/Pk${mode}samples*`
    ${FSLDIR}/bin/fslmaths ${out_dir}/merged_Pk${mode}samples -Tmean ${out_dir}/mean_Pk${mode}samples
    ${FSLDIR}/bin/make_dyadic_vectors ${out_dir}/merged_Pth${mode}samples ${out_dir}/merged_Pph${mode}samples ${out_dir}/mean_Pk${mode}samples ${out_dir}/Pdyads${mode}
    mode=$(($mode + 1))
done


echo Removing intermediate files

if [ `${FSLDIR}/bin/imtest ${out_dir}/merged_th1samples` -eq 1 ];then
  if [ `${FSLDIR}/bin/imtest ${out_dir}/merged_ph1samples` -eq 1 ];then
    if [ `${FSLDIR}/bin/imtest ${out_dir}/merged_f1samples` -eq 1 ];then
      rm -rf ${out_dir}/diff_slices
      rm -f ${out_dir}/dataLR_slice_*
      rm -f ${out_dir}/dataHR_newslice_* 
       rm -f ${out_dir}/nodif_brain_maskLR_slice_*
      if [ `${FSLDIR}/bin/imtest ${out_dir}/grad_devLR_slice_0000` -eq 1 ];then
	  rm -f ${out_dir}/grad_devLR_slice_*
      fi	  
      if [ `${FSLDIR}/bin/imtest ${out_dir}/grad_devHR_newslice_0000` -eq 1 ];then
	  rm -f ${out_dir}/grad_devHR_newslice_*
      fi	  
    fi
  fi
fi

echo Done
