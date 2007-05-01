
#   FSL interface for FDT (BEDPOST and ProbTrack)
#
#   Timothy Behrens, Heidi Johansen-Berg, Dave Flitney and Matthew Webster FMRIB Image Analysis Group
#
#   Copyright (C) 2006 University of Oxford
#
#   TCLCOPYRIGHT

#TO DO replace -filetypes * with -filetypes { } for directory selectors
source [ file dirname [ info script ] ]/fslstart.tcl
option add *FileEntry*Entry*width 35
set TCLPATH [file dirname [ info script ] ]
regsub tcl $TCLPATH bin BINPATH
regsub tcl $TCLPATH doc/fdt HTMLPATH

set VERSION "1.1x"

proc mm_to_voxels { X Y Z mask } {

    global FSLDIR

    upvar $X cX
    upvar $Y cY
    upvar $Z cZ


    set vcX [ exec sh -c "echo $cX $cY $cZ | $FSLDIR/bin/std2imgcoord -img $mask -vox - | awk '{print \$1}'" ]    
    set vcY [ exec sh -c "echo $cX $cY $cZ | $FSLDIR/bin/std2imgcoord -img $mask -vox - | awk '{print \$2}'" ] 
    set vcZ [ exec sh -c "echo $cX $cY $cZ | $FSLDIR/bin/std2imgcoord -img $mask -vox - | awk '{print \$3}'" ] 	
    set cX $vcX
    set cY $vcY
    set cZ $vcZ
}

proc fdt:dialog { w tclstartupfile } {

    global eddy bedpost registration dtifit probtrack HTMLPATH FSLDIR VERSION INMEDX VARS
    set probtrack(tool) "probtrack"

    if [winfo exists $w] {
        wm deiconify $w
        raise $w
        return
    }

    toplevel $w
    wm title $w "FDT - FMRIB's Diffusion Toolbox $VERSION"
    wm iconname $w "FDT"
    wm iconbitmap $w @$FSLDIR/tcl/fmrib.xbm

    #-------- Stage and Mode Options -------- 

    frame $w.tool
    optionMenu2 $w.tool.menu probtrack(tool) -command "fdt:select_tool $w" eddy_current "Eddy current correction" bedpost "BEDPOSTX Estimation of diffusion parameters"  registration "Registration" probtrack "ProbTrackX Probabilistic tracking" xutilssx "----------------------------------------------------" dtifit "DTIFit Reconstruct diffusion tensors" 
    $w.tool.menu.menu entryconfigure 4 -state disabled -background black
    pack $w.tool.menu -side left -pady 3 -padx 6 -anchor nw

    #-------- Tool Options... -------- 

    frame $w.opts

    #------- Registration --------

    frame $w.registration

    proc registration_set_directory { w dirname } {
	global registration

	set struct [ file join $dirname struct_brain ]

	if { [ imtest $struct ] } {
	    set registration(struct_image) $struct
	} else {
	    set registration(struct_image) ""
	}
    }

    FileEntry $w.registration.directory -textvariable registration(directory) -label "BEDPOST directory:" -title "Choose directory" -filetypes * -command "registration_set_directory $w" 

    frame       $w.registration.struct
    checkbutton $w.registration.struct.yn -variable registration(struct_yn) -command "registration_packframe $w"
    label       $w.registration.struct.lb -text "Main structural image"
    TitleFrame  $w.registration.struct.tf -text "Main structural image" 
    optionMenu2 $w.registration.struct.tf.search registration(struct_search) 0 "No search" 90 "Normal search" 180 "Full search"
    optionMenu2 $w.registration.struct.tf.dof registration(struct_dof)   6 "6 DOF" 7 "7 DOF" 9 "9 DOF" 12 "12 DOF"  
    optionMenu2 $w.registration.struct.tf.costfn registration(struct_costfn) corratio "Correlation ratio" mutualinfo "Mutual information"
    FileEntry   $w.registration.struct.tf.file -textvariable registration(struct_image) -filetypes IMAGE -width 45
    pack $w.registration.struct.tf.file -side top -in [ $w.registration.struct.tf getframe ]
    pack $w.registration.struct.tf.search $w.registration.struct.tf.dof $w.registration.struct.tf.costfn -side left  -in [ $w.registration.struct.tf getframe ]
    set registration(struct_costfn) mutualinfo
    set registration(struct_dof) 12
    set registration(struct_search) 90
    set registration(struct_yn) 0

    frame       $w.registration.standard
    checkbutton $w.registration.standard.yn -variable registration(standard_yn)  -command "registration_packframe $w"
    TitleFrame  $w.registration.standard.tf -text "Standard space"
    label       $w.registration.standard.lb -text "Standard space"
    optionMenu2 $w.registration.standard.tf.search registration(standard_search) 0 "No search" 90 "Normal search" 180 "Full search"
    optionMenu2 $w.registration.standard.tf.dof registration(standard_dof)   6 "6 DOF" 7 "7 DOF" 9 "9 DOF" 12 "12 DOF"  
    optionMenu2 $w.registration.standard.tf.costfn registration(standard_costfn) corratio "Correlation ratio" mutualinfo "Mutual information"
    FileEntry   $w.registration.standard.tf.file -textvariable registration(standard_image) -filetypes IMAGE -width 45 
    pack $w.registration.standard.tf.file -side top -in [ $w.registration.standard.tf getframe ]
    pack $w.registration.standard.tf.search $w.registration.standard.tf.dof $w.registration.standard.tf.costfn -side left -in [ $w.registration.standard.tf getframe ]
    set registration(standard_yn) 1
    set registration(standard_dof) 12
    set registration(standard_search) 90
    set registration(standard_image) [ file join ${FSLDIR} etc standard avg152T1_brain ]

    pack $w.registration.directory $w.registration.struct $w.registration.standard -side top -padx 3 -pady 3 -anchor w

    proc registration_packframe { w } {
       global registration
       pack forget $w.registration.struct.yn $w.registration.struct.tf $w.registration.struct.yn $w.registration.struct.lb
       pack forget $w.registration.standard.yn $w.registration.standard.tf $w.registration.standard.yn $w.registration.standard.lb
       if {$registration(struct_yn)} { pack $w.registration.struct.yn $w.registration.struct.tf -side left -anchor w } else { pack $w.registration.struct.yn  $w.registration.struct.lb -side left -anchor w}
       if {$registration(standard_yn)} { pack $w.registration.standard.yn $w.registration.standard.tf -side left -anchor w } else { pack $w.registration.standard.yn  $w.registration.standard.lb -side left -anchor w}
    }
    
    registration_packframe $w
    #------- ECC --------
    frame $w.ecc

    proc ecc_update_files { w filename } {
	global eddy
	set eddy(output) [ file join [file dirname $eddy(input)] data ]
    }    

    FileEntry $w.ecc.input -textvariable eddy(input) -label "Diffusion weighted data:" -title "Choose diffusion weighted image" -filetypes IMAGE -command "ecc_update_files $w"

    FileEntry $w.ecc.output -textvariable eddy(output) 	-label "Corrected output data:" -title  "Choose output image name" -filetypes IMAGE -command "ecc_update_files $w"

   set eddy(refnum) 0
   LabelSpinBox  $w.ecc.refnum -label "Reference volume"  -textvariable eddy(refnum) -range {0 100 1 } -width 6 

    pack $w.ecc.input $w.ecc.output $w.ecc.refnum -side top -padx 3 -pady 3 -expand yes -anchor w

   #------- DTIFit --------

    frame $w.dtifit

    FileEntry $w.dtifit.directory -textvariable dtifit(directory) -label  "Input directory:" -title "Choose directory" -command "set_working_directory dtifit(cwd)"

    proc dtifit_toggle_expert { w } {
	global dtifit

	if { $dtifit(expert_yn) } {
	    pack forget $w.dtifit.directory
	    pack $w.dtifit.expert -in $w.dtifit -after $w.dtifit.expert_yn
	} else {
	    pack forget $w.dtifit.expert
	    pack $w.dtifit.directory -in $w.dtifit -before $w.dtifit.expert_yn
	}
    }

    checkbutton $w.dtifit.expert_yn -text "Specify input files manually" \
	-variable dtifit(expert_yn) -command "dtifit_toggle_expert $w"

    frame $w.dtifit.expert

    proc set_working_directory { cwd filename } {
	global dtifit
	set dirname [file dirname $filename]
	puts "switching from $dtifit(cwd) to $dirname" 
	set dtifit(cwd) $dirname
    }

    proc dtifit_update_files { w filename } {
	global dtifit

	set dtifit(output) [ file join [file dirname $dtifit(input)] dti ]
	set_working_directory dtifit(cwd) $dtifit(input)
    }
    
    set dtifit(cwd) [ pwd ]

#All the below orignally had -directory $dtifit(cwd) 
    option add *dtifit.expert.FileEntry*labf*width 27
    FileEntry $w.dtifit.expert.input -textvariable dtifit(input) -label  "Diffusion weighted data:" -title "Choose diffusion weighted image" -filetypes IMAGE -command "dtifit_update_files $w" 
    FileEntry $w.dtifit.expert.mask -textvariable dtifit(mask) -label "BET binary brain mask:" -title "Choose BET brain mask file" -filetypes IMAGE -command "set_working_directory dtifit(cwd)"
    FileEntry $w.dtifit.expert.output -textvariable dtifit(output) -label "Output basename:" -title  "Choose output base name" -filetypes * -command "set_working_directory dtifit(cwd)"
    FileEntry $w.dtifit.expert.bvecs -textvariable dtifit(bvecs) -label "Gradient directions:" -title  "Choose bvecs file" -filetypes * -command "set_working_directory dtifit(cwd)"
    FileEntry $w.dtifit.expert.bvals -textvariable dtifit(bvals) -label  "b values:" -title  "Choose bvals file" -command "set_working_directory dtifit(cwd)"

    pack $w.dtifit.expert.input $w.dtifit.expert.mask $w.dtifit.expert.output \
	$w.dtifit.expert.bvecs $w.dtifit.expert.bvals \
	-side top -padx 3 -pady 3 -expand yes -anchor w

    pack $w.dtifit.directory $w.dtifit.expert_yn -side top -padx 3 -pady 3 -expand yes -anchor w

    #------- BEDPOST --------

    frame $w.bedpost

    FileEntry $w.bedpost.directory -textvariable bedpost(directory) -label "Input directory:" -title "Choose directory" -filetypes * -command "set_working_directory dtifit(cwd)"

   collapsible frame $w.bedpost.advanced -title "Advanced Options"
   set bedpost(nfibres) 2
   set bedpost(weight)  1
   set bedpost(burnin)  1000
   LabelSpinBox  $w.bedpost.advanced.nfibres -label "Fibres  "  -textvariable bedpost(nfibres) -range {1 1000000000 1 } 
   LabelSpinBox  $w.bedpost.advanced.weight  -label "Weight "  -textvariable bedpost(weight) -range {0.0 100000000.0 1 } 
   LabelSpinBox  $w.bedpost.advanced.burnin  -label "Burn In"  -textvariable bedpost(burnin) -range {1 1000000000 1 } 

    set bedpost(ecc_yn) 0
    pack $w.bedpost.advanced.nfibres $w.bedpost.advanced.weight $w.bedpost.advanced.burnin -in $w.bedpost.advanced.b -anchor w
    pack $w.bedpost.directory $w.bedpost.advanced -side top -padx 3 -pady 3 -expand yes -anchor w

    #-------- ProbTrack -------- 
    NoteBook $w.probtrack -bd 2 -tabpady {5 10} -arcradius 3
    $w.probtrack insert 0 data -text "Data"
    $w.probtrack insert 1 options -text "Options"
    #-------- Mode specific option --------
    frame $w.data
    FileEntry $w.data.directory -textvariable probtrack(bedpost_dir) -label "BEDPOST directory" -title "Choose BEDPOST directory" -filetypes * -command "probtrack_update_files $w"

    TitleFrame  $w.data.seed -text "Seed Space"
    optionMenu2 $w.data.seed.menu probtrack(mode) -command "fdt:probtrack_mode $w" simple "Single voxel" seedmask "Single mask" network "Multiple masks"
    set probtrack(x) 0
    set probtrack(y) 0
    set probtrack(z) 0
    set probtrack(units) vox
    #Co-ordinate edit frame
    frame $w.data.seed.voxel
    LabelSpinBox $w.data.seed.voxel.x -label "X" -textvariable probtrack(x) -range {-1000000 1000000 1 } 
    LabelSpinBox $w.data.seed.voxel.y -label "Y" -textvariable probtrack(y) -range {-1000000 1000000 1 } 
    LabelSpinBox $w.data.seed.voxel.z -label "Z" -textvariable probtrack(z) -range {-1000000 1000000 1 } 
    radiobutton $w.data.seed.voxel.vox -text "vox" -value vox -variable probtrack(units)
    radiobutton $w.data.seed.voxel.mm  -text "mm"  -value mm  -variable probtrack(units)
    FileEntry $w.data.seed.reference -textvariable probtrack(reference) -label "Seed reference image:" -title "Choose reference image" -filetypes IMAGE 

    option add *seed*FileEntry*labf*width 24

    frame  $w.data.seed.ssf
    checkbutton $w.data.seed.ssf.ssd -text "Seed space is not diffusion" -variable probtrack(usereference_yn)  -command " pack forget $w.data.seed.ssf.xfm  ; if { \$probtrack(usereference_yn) } { pack $w.data.seed.ssf.xfm } ; $w.probtrack compute_size"
    FileEntry $w.data.seed.ssf.xfm -textvariable probtrack(xfm)  -label "Select Seed to diff transform" -title "Select seed-space to DTI-space transformation matrix" -filetypes *
    pack $w.data.seed.ssf.ssd -side top -anchor nw

    


    frame  $w.data.seed.bcf
    checkbutton $w.data.seed.bcf.bc -text "Blind Classification:" -variable probtrack(bcyn) 
    pack $w.data.seed.bcf.bc -side top -anchor w

    pack $w.data.seed.voxel.x $w.data.seed.voxel.y $w.data.seed.voxel.z $w.data.seed.voxel.vox $w.data.seed.voxel.mm -side left -padx 2
    pack $w.data.seed.voxel $w.data.seed.ssf -in $w.data.seed.f -side bottom -anchor w -pady 2
    pack $w.data.seed.menu $w.data.seed.reference -in $w.data.seed.f -side left -anchor w -pady 2

    TitleFrame $w.data.seed.target -text "Multiple Masks"    
    listbox $w.data.seed.targets -height 6 -width 50 -yscrollcommand "$w.data.seed.sb set"
    scrollbar $w.data.seed.sb -command "$w.data.seed.targets yview " 
    frame $w.data.seed.tb
    button $w.data.seed.tb.add -text "Add Image"  -command "feat_file:setup_dialog $w a a a [namespace current] IMAGE {Select File} {fdt_add $w $w.data.seed.targets} {}" 
    button $w.data.seed.tb.del -text "Remove Image"  -command "fdt_sub $w $w.data.seed.targets" 
    button $w.data.seed.tb.imp -text "Load List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_imp $w $w.data.seed.targets} {}"
    button $w.data.seed.tb.exp -text "Save List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_exp $w $w.data.seed.targets} {}"
    pack $w.data.seed.tb.add $w.data.seed.tb.del $w.data.seed.tb.imp $w.data.seed.tb.exp -side left
    pack $w.data.seed.tb -in [$w.data.seed.target getframe ] -side bottom  -expand yes -fill x -anchor w -padx 3 -pady 3
    pack $w.data.seed.targets $w.data.seed.sb -in [$w.data.seed.target getframe ] -side left  -expand yes -fill y -anchor w -padx 3 -pady 3
    
    TitleFrame  $w.data.targets -text "Targets"

    proc fdt_add { w listbox filename } {
    set filename [ fix_cygwin_filename $filename ]
    $listbox insert end $filename
    }

    proc fdt_sub { w listbox} {
    set count 0
    foreach file [ $listbox get 0 end ] {
	if { [ $listbox selection includes $count ] == 1 } {
	    $listbox delete $count
	    incr count -1
	}
	incr count
    } 
    }

    proc fdt_imp { w listbox filename } {
    if { ![ file exists $filename ] } {
	MxPause "Warning: Bad or missing file!"
	return
    }
    set fd [ open $filename ]
    $listbox  delete 0 end
    while { [ gets $fd file ] >= 0 } {
	$listbox insert end $file
    }
    close $fd
    }

    proc fdt_exp { w listbox filename } {
    set fd [ open $filename w ]
    foreach file [ $listbox get 0 end ] {
	puts $fd $file
    }
    close $fd
    }

    frame $w.data.targets.wf    
    checkbutton $w.data.targets.wf.sct -text "Set Waypoints" -variable probtrack(waypoint_yn)  -command " pack forget $w.data.targets.wf.tf ; if { \$probtrack(waypoint_yn) } { pack $w.data.targets.wf.tf } ; $w.probtrack compute_size"
    TitleFrame $w.data.targets.wf.tf -text "Waypoints"    
    listbox $w.data.targets.wf.tf.targets -height 6 -width 50 -yscrollcommand "$w.data.targets.wf.tf.sb set"
    scrollbar $w.data.targets.wf.tf.sb -command "$w.data.targets.wf.tf.targets yview " 
    frame $w.data.targets.wf.tf.tb
    button $w.data.targets.wf.tf.tb.add -text "Add Image"  -command "feat_file:setup_dialog $w a a a [namespace current] IMAGE {Select File} {fdt_add $w $w.data.targets.wf.tf.targets} {}"
    button $w.data.targets.wf.tf.tb.del -text "Remove Image"  -command "fdt_sub $w $w.data.targets.wf.tf.targets" 
    button $w.data.targets.wf.tf.tb.imp -text "Load List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_imp $w $w.data.targets.wf.tf.targets} {}"
    button $w.data.targets.wf.tf.tb.exp -text "Save List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_exp $w $w.data.targets.wf.tf.targets} {}"
    pack $w.data.targets.wf.tf.tb.add $w.data.targets.wf.tf.tb.del $w.data.targets.wf.tf.tb.imp $w.data.targets.wf.tf.tb.exp -side left
    pack $w.data.targets.wf.tf.tb -in [$w.data.targets.wf.tf getframe ] -side bottom  -expand yes -fill x -anchor w -padx 3 -pady 3
    pack $w.data.targets.wf.tf.targets $w.data.targets.wf.tf.sb -in [$w.data.targets.wf.tf getframe ] -side left  -expand yes -fill y -anchor w -padx 3 -pady 3
    pack  $w.data.targets.wf.sct -side top -anchor nw
    pack $w.data.targets.wf 
  
    option add *targets*Checkbutton*width 18
    option add *targets*Checkbutton*anchor w
    frame  $w.data.targets.ef
    checkbutton $w.data.targets.ef.srt -text "Set Exclusion targets" -variable probtrack(exclude_yn)  -command " pack forget $w.data.targets.ef.rubbish ; if { \$probtrack(exclude_yn) } { pack $w.data.targets.ef.rubbish } ; $w.probtrack compute_size"
    FileEntry $w.data.targets.ef.rubbish -textvariable probtrack(exclude) -title "Select exclusion image" -filetypes IMAGE
    pack $w.data.targets.ef.srt -side left

    frame  $w.data.targets.sf
    checkbutton $w.data.targets.sf.sst -text "Set Termination targets" -variable probtrack(terminate_yn)  -command " pack forget $w.data.targets.sf.stop ; if { \$probtrack(terminate_yn) } { pack $w.data.targets.sf.stop } ; $w.probtrack compute_size"
    FileEntry $w.data.targets.sf.stop -textvariable probtrack(stop) -title "Select termination image" -filetypes IMAGE
    pack $w.data.targets.sf.sst -side left

    frame $w.data.targets.cf    
    checkbutton $w.data.targets.cf.sct -text "Set Classification targets" -variable probtrack(classify_yn)  -command " pack forget $w.data.targets.cf.tf ; if { \$probtrack(classify_yn) } { pack $w.data.targets.cf.tf } ; $w.probtrack compute_size"
    TitleFrame $w.data.targets.cf.tf -text "Classification"    
    listbox $w.data.targets.cf.tf.targets -height 6 -width 50 -yscrollcommand "$w.data.targets.cf.tf.sb set"
    scrollbar $w.data.targets.cf.tf.sb -command "$w.data.targets.cf.tf.targets yview " 
    frame $w.data.targets.cf.tf.tb
    button $w.data.targets.cf.tf.tb.add -text "Add Image"  -command "feat_file:setup_dialog $w a a a [namespace current] IMAGE {Select File} {fdt_add $w $w.data.targets.cf.tf.targets} {}"
    button $w.data.targets.cf.tf.tb.del -text "Remove Image"  -command "fdt_sub $w $w.data.targets.cf.tf.targets" 
    button $w.data.targets.cf.tf.tb.imp -text "Load List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_imp $w $w.data.targets.cf.tf.targets} {}"
    button $w.data.targets.cf.tf.tb.exp -text "Save List" -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_exp $w $w.data.targets.cf.tf.targets} {}"
    pack $w.data.targets.cf.tf.tb.add $w.data.targets.cf.tf.tb.del $w.data.targets.cf.tf.tb.imp $w.data.targets.cf.tf.tb.exp -side left
    pack $w.data.targets.cf.tf.tb -in [$w.data.targets.cf.tf getframe ] -side bottom  -expand yes -fill x -anchor w -padx 3 -pady 3
    pack $w.data.targets.cf.tf.targets $w.data.targets.cf.tf.sb -in [$w.data.targets.cf.tf getframe ] -side left  -expand yes -fill y -anchor w -padx 3 -pady 3
    pack  $w.data.targets.cf.sct -side top -anchor nw
    pack $w.data.targets.cf 


    pack $w.data.targets.wf $w.data.targets.ef $w.data.targets.sf $w.data.targets.cf -in $w.data.targets.f -anchor w




    set probtrack(xfm) ""
    set probtrack(basename) "merged"
    set probtrack(mask) "nodif_brain_mask"

    proc probtrack_update_files { w filename } {
	global probtrack
	global FSLDIR

	if { ($probtrack(bedpost_dir) != "") && ($probtrack(reference) != "") } {
	    set probtrack(output) \
		[ file join $probtrack(bedpost_dir) [ file tail [ exec $FSLDIR/bin/remove_ext $probtrack(reference) ] ] ]
	}
    }

    FileEntry $w.data.dir -textvariable probtrack(output) -label  "Output directory:" -title  "Name the output directory" -filetypes *

    pack $w.data.directory $w.data.seed $w.data.targets $w.data.dir -padx 3 -pady 3 -anchor nw

    pack $w.data -in  [$w.probtrack getframe data] -padx 3 -pady 3 -anchor nw -expand yes -fill both

    #-------- ...Options... --------
    TitleFrame $w.options -text "Basic Options"

    checkbutton $w.options.verbose -text "Verbose" -variable probtrack(verbose_yn)
    
    set probtrack(nparticles) 5000
    LabelSpinBox $w.options.nparticles -label  "Number of samples" -textvariable probtrack(nparticles) -range {1 1e24 100 } -width 6 

    set probtrack(curvature) 0.2
    LabelSpinBox $w.options.curvature -label "Curvature threshold" -textvariable probtrack(curvature) -range {0.0 1.0 0.01 }

    set probtrack(loopcheck_yn) 1
    checkbutton $w.options.loopcheck -text "Loopcheck" -variable probtrack(loopcheck_yn)

    collapsible frame $w.advanced -title "Advanced Options" -command "$w.probtrack compute_size; set dummy"

    set probtrack(nsteps) 2000
    LabelSpinBox $w.advanced.nsteps -label "Maximum number of steps" -textvariable probtrack(nsteps) -range {2 1000000 10 } -width 6

    set probtrack(steplength) 0.5
    LabelSpinBox $w.advanced.steplength -label "Step length (mm)" -textvariable probtrack(steplength) -range {0 10000 0.1} 

    set probtrack(modeuler_yn) 0
    checkbutton $w.advanced.modeuler -text "Use modified Euler streamlining" -variable probtrack(modeuler_yn)

    set probtrack(pd) 0
    checkbutton $w.advanced.pd -text "Use Distance correction" -variable probtrack(pd)

    set probtrack(usef_yn) 0
    checkbutton $w.advanced.usef -text "Use anisotropy to constrain tracking" -variable probtrack(usef_yn)

    pack $w.advanced.modeuler $w.advanced.nsteps $w.advanced.steplength $w.advanced.usef $w.advanced.pd -in $w.advanced.b  -side top -pady 3 -padx 6 -anchor nw

    pack \
	$w.options.nparticles \
	$w.options.curvature \
	$w.options.verbose \
	$w.options.loopcheck \
	-in [$w.options getframe ] -side top -pady 3 -padx 6 -anchor nw

    pack $w.options $w.advanced -in [$w.probtrack getframe options] -side top -pady 3 -padx 6 -anchor nw -expand yes -fill both

    #-------- Buttons --------

    frame $w.btns
    frame $w.btns.b -relief raised -borderwidth 1
    
    button $w.apply     -command "fdt:apply $w keep" \
        -text "Go" -width 5
    bind $w.apply <Return> {
        [winfo toplevel %W].apply invoke}
 
    button $w.cancel    -command "fdt:destroy $w" \
        -text "Exit" -width 5
    bind $w.cancel <Return> {
        [winfo toplevel %W].cancel invoke}

    button $w.help -command "FmribWebHelp file: $HTMLPATH/index.html" \
	    -text "Help" -width 5
    bind $w.help <Return> {
	[winfo toplevel %W].help invoke}
 
    pack $w.btns.b -side bottom -fill x -padx 3 -pady 5
    pack $w.apply $w.cancel $w.help -in $w.btns.b \
        -side left -expand yes -padx 3 -pady 10 -fill y
 
    pack $w.tool $w.opts $w.btns -side top -expand yes -fill both
    
 

    $w.probtrack raise data 
    fdt:select_tool $w 
    set probtrack(mode) simple
    fdt:probtrack_mode $w

    update idletasks
    if { $tclstartupfile != "" } {
	puts "Reading $tclstartupfile"
	source $tclstartupfile
	fdt:select_tool $w 
	fdt:probtrack_mode $w
    }
}

proc fdt:probtrack_mode { w } {
    global probtrack

    pack forget $w.data.seed.voxel $w.data.seed.ssf $w.data.seed.menu $w.data.seed.reference $w.data.seed.bcf $w.data.seed.target $w.data.targets.cf
    $w.data.dir configure -label  "Output directory:" -title  "Name the output directory" -filetypes *
    switch -- $probtrack(mode) {
  	simple {
                     pack $w.data.seed.ssf $w.data.seed.voxel -in $w.data.seed.f -side bottom -anchor w -pady 2
                     pack $w.data.seed.menu $w.data.seed.reference -in $w.data.seed.f -side left -anchor w -pady 2
                     $w.data.seed.reference configure -label "Seed reference image:" -title "Choose reference image" 
                     $w.data.dir configure -label  "Output file:" -title  "Name the output file" -filetypes IMAGE
    	}
	seedmask {
                     pack $w.data.seed.ssf $w.data.seed.bcf -in $w.data.seed.f -side bottom -anchor w -pady 2
                     pack $w.data.seed.menu $w.data.seed.reference -in $w.data.seed.f -side left -anchor w -pady 2
                     pack $w.data.targets.cf -in $w.data.targets.f -anchor w
                     $w.data.seed.reference configure -label "Mask image:" -title "Choose mask image" 
  	}
	network {
                     pack  $w.data.seed.target $w.data.seed.ssf $w.data.seed.menu -in $w.data.seed.f -side bottom -anchor w -pady 2
	}
    }
    $w.probtrack compute_size
}

proc fdt:select_tool { w } {
    global probtrack
    pack forget $w.ecc
    pack forget $w.probtrack
    pack forget $w.bedpost
    pack forget $w.registration
    pack forget $w.dtifit
    if {$probtrack(tool) == "bedpost"} { 
	pack $w.bedpost -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$probtrack(tool) == "probtrack"}  { 
	pack $w.probtrack -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$probtrack(tool) == "dtifit"}  { 
	pack $w.dtifit -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$probtrack(tool) == "eddy_current"}  { 
	pack $w.ecc -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    } elseif {$probtrack(tool) == "registration"} {
	pack $w.registration -in $w.opts -side top -padx 3 -pady 3 -anchor nw
    }
}
proc fdt_monitor_short { w cmd } {
    global debugging OSFLAVOUR

    puts "$cmd"

    if { $OSFLAVOUR != "cygwin" } {
	set oldcursor [ $w configure -cursor { watch red white } ]

	catch {
	    update idletasks
	    if { ! $debugging } {
		set fd [ open "|$cmd" r ]
#		set fd [ open "|qrsh -V -now n -q short.q $cmd" r ]
		while { ( [ gets $fd line ] >= 0 ) } {
		    update idletasks
		    puts $line
		}
		close $fd
	    }
	} junk

	$w configure -cursor $oldcursor

    } else {
	catch { exec sh -c $cmd } junk
    }

    if { $junk != "" } {
	MxPause "Errors: $junk"
    } 

    puts "Done!"
}

proc fdt_monitor { w cmd } {
    global debugging OSFLAVOUR

    puts "$cmd"

    if { $OSFLAVOUR != "cygwin" } {
	set oldcursor [ $w configure -cursor { watch red white } ]

	catch {
	    update idletasks
	    if { ! $debugging } {
		set fd [ open "|$cmd" r ]
#		set fd [ open "|qrsh -V -now n -q long.q $cmd" r ]
		while { ( [ gets $fd line ] >= 0 ) } {
		    update idletasks
		    puts $line
		}
		close $fd
	    }
	} junk

	$w configure -cursor $oldcursor

    } else {
	catch { exec sh -c $cmd } junk
    }

    if { $junk != "" } {
	MxPause "Errors: $junk"
    } 

    puts "Done!"
}

proc fdt:apply { w dialog } {

    global probtrack
    global BINPATH
    global FSLDIR

    switch -- $probtrack(tool) {
	eddy_current {
	    global eddy

	    set errorStr ""
	    if { $eddy(input) == "" } { set errorStr "You need to specify the input image! " }
	    if { $eddy(output) == "" } { set errorStr "$errorStr You need to specify an output image!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    #	    check output!=input
	    set canwrite 1
	    if { $eddy(input) == $eddy(output) } {
		set canwrite [ YesNoWidget "Output and input images have the same name. Overwrite input?" Yes No ]
	    }
	    if { $canwrite } {
		fdt_monitor $w "${FSLDIR}/bin/eddy_correct $eddy(input) $eddy(output) $eddy(refnum)"
	    }
	}
	dtifit {
	    global dtifit

	    if { ! $dtifit(expert_yn) } {
		set dtifit(input)  [ file join $dtifit(directory) data ]
		set dtifit(output) [ file join $dtifit(directory) dti ]
		set dtifit(mask)   [ file join $dtifit(directory) nodif_brain_mask ]
		set dtifit(bvecs)  [ file join $dtifit(directory) bvecs ]
		set dtifit(bvals)  [ file join $dtifit(directory) bvals ]
	    }

	    set errorStr ""
	    if { $dtifit(directory) == "" && ! $dtifit(expert_yn) } { set errorStr "You must specify the input directory!" }
	    if { $dtifit(input) == "" } { set errorStr "You need to specify the diffusion weighted data image!" }
	    if { $dtifit(output) == "" } { set errorStr "$errorStr You need to specify the output basename!" }
	    if { $dtifit(mask) == "" } { set errorStr "$errorStr You need to specify a mask image!" }
	    if { $dtifit(bvecs) == "" } { set errorStr "$errorStr Please select a gradient directions file!" }
	    if { $dtifit(bvals) == "" } { set errorStr "$errorStr Please select a b values file!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    set canwrite 1
	    if { [file exists $dtifit(output) ] } {
		set canwrite [ YesNoWidget "Overwrite $dtifit(output)?" Yes No ]
	    }
	    if { $canwrite } {
		fdt_monitor_short $w "${FSLDIR}/bin/dtifit --data=$dtifit(input) --out=$dtifit(output) --mask=$dtifit(mask) --bvecs=$dtifit(bvecs) --bvals=$dtifit(bvals)"
	    }
	}
	bedpost {
	    global bedpost

	    set errorStr ""
	    if { $bedpost(directory) == ""  } { set errorStr "You must specify the bedpost directory!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    set canwrite 1
	    if { [file exists ${bedpost(directory)}.bedpost ] } {
		set canwrite [ YesNoWidget "Overwrite ${bedpost(directory)}.bedpost?" Yes No ]
		if { $canwrite } {
		    puts "rm -rf ${bedpost(directory)}.bedpost"
		    catch { exec rm -rf ${bedpost(directory)}.bedpost } errmsg
		}
	    }
	    if { $canwrite } {
		puts "bedpostX $bedpost(directory) -n $bedpost(nfibres) -w $bedpost(weight)  -b $bedpost(burnin)"
		fdt_monitor $w "${FSLDIR}/bin/bedpostX $bedpost(directory) -n $bedpost(nfibres) -w $bedpost(weight)  -b $bedpost(burnin)"
	    }
	}
	probtrack {
	    global probtrack
	    set errorStr ""
	    if { $probtrack(bedpost_dir) == ""  } { set errorStr "You must specify the bedpost directory!" }
	    if { $probtrack(mode) != "network" && $probtrack(reference) == "" } { set errorStr "$errorStr You must specify a reference image" } 
	    if { $probtrack(exclude_yn) && $probtrack(exclude) == "" } { set errorStr "$errorStr You must specify the exclusion mask!" }
            if { $probtrack(terminate_yn) && $probtrack(stop) == ""} { set errorStr "$errorStr You must specify the termination mask!" }
	    if { $probtrack(output) == ""  } { set errorStr "$errorStr You must specify the output basename!" }
	    set flags ""
	    if { $probtrack(verbose_yn) == 1 } { set flags "$flags -V 1" }
	    if { $probtrack(loopcheck_yn) == 1 } { set flags "$flags -l" }
	    if { $probtrack(usef_yn) == 1 } { set flags "$flags -f" }
	    if { $probtrack(modeuler_yn) == 1 } { set flags "$flags --modeuler" }
	    set flags "$flags -c $probtrack(curvature) -S $probtrack(nsteps) --steplength=$probtrack(steplength) -P $probtrack(nparticles)"
             
            if { $probtrack(pd) } { set flags "$flags --pd"  }

	    set tn [open "| $BINPATH/tmpnam"]
	    gets $tn filebase
	    close $tn
	    set logfile "${filebase}_log.tcl"
	    set log [open "$logfile" w]
	    puts $log "set tool $probtrack(tool)"
	    set copylog ""

	    if { $probtrack(usereference_yn) } {
		set flags "$flags --xfm=$probtrack(xfm)"
      		puts $log "set probtrack(xfm) $probtrack(xfm)"
	    }

	    if { $probtrack(exclude_yn) == 1 } {
		set flags "$flags --avoid=$probtrack(exclude)"
		puts $log "set probtrack(exclude) $probtrack(exclude)"
	    }

	    if { $probtrack(terminate_yn) == 1 } {
		set flags "$flags --stop=$probtrack(stop)"
		puts $log "set probtrack(stop) $probtrack(stop)"
	    }

	    if { $errorStr != "" } {
       		MxPause $errorStr
       		return
      	    }
	    set canwrite 1
      	    if { [ file exists $probtrack(output) ] } {
      		set canwrite [  YesNoWidget "Overwrite $probtrack(output)?" Yes No ]
	    }
       	    if { $canwrite } {
       		puts "rm -rf $probtrack(output)"
       		exec rm -rf $probtrack(output)
	        puts "mkdir -p $probtrack(output)"
		exec mkdir -p $probtrack(output)
       	    }

	    set flags "$flags --forcedir --opd -s $probtrack(bedpost_dir)/merged -m $probtrack(bedpost_dir)/nodif_brain_mask  --dir=$probtrack(output)" 
    	    foreach entry {bedpost_dir xfm mode exclude_yn usereference_yn verbose_yn loopcheck_yn modeuler_yn curvature nsteps steplength nparticles} {
		puts $log "set probtrack($entry) $probtrack($entry)"
	    }
            switch $probtrack(mode) {
	       simple { 
		    set fd [ open "${filebase}_coordinates.txt" w ]
		    set x $probtrack(x)
		    set y $probtrack(y)
		    set z $probtrack(z)
		    if { $probtrack(units) == "mm" } {
			if { $probtrack(reference) != "" } {
			    mm_to_voxels x y z $probtrack(reference)
			} else {
			    mm_to_voxels x y z [ file join $probtrack(bedpost_dir) nodif_brain_mask ]
			}			    
			puts $fd "$x $y $z"
			puts "$probtrack(x) $probtrack(y) $probtrack(z) (mm) -> $x $y $z (voxels)"
		    } else {
			puts $fd "$probtrack(x) $probtrack(y) $probtrack(z)"
		    }
		    close $fd
 		    puts $log "set probtrack(x) $probtrack(x)"
		    puts $log "set probtrack(y) $probtrack(y)"
		    puts $log "set probtrack(z) $probtrack(z)"
		    puts $log "set probtrack(units) $probtrack(units)"
                    set flags "--mode=single --seedref=$probtrack(reference) -o $probtrack(output) -x ${filebase}_coordinates.txt $flags"
	       } 
               seedmask {
		     if { $probtrack(bcyn) } { 
                       fdt_monitor_short $w "${FSLDIR}/bin/flirt -in $probtrack(bedpost_dir)/seed_brain_mask -ref $probtrack(bedpost_dir)/seed_brain_mask -applyisoxfm 5 -out $probtrack(output)/lowresmask"
                       fdt_monitor_short $w "${FSLDIR}/bin/avwmaths  $probtrack(output)/lowresmask -thr 0.5 -bin  $probtrack(output)/lowresmask"
                       set flags "$flags --lrmask=$probtrack(output)/lowresmask --omatrix2" 
                   }
                   set flags "--mode=seedmask -x $probtrack(reference) $flags"  
	       }
	       network {
                   fdt_exp w $w.data.seed.targets $probtrack(output)/masks.txt
		   set flags "--network --mode=seedmask -x $probtrack(output)/masks.txt $flags"
	       }
	    }

       	    if { $canwrite } {
       		set copylog "fdt.log"

	        if { $probtrack(waypoint_yn) == 1 } {
                    fdt_exp w $w.data.targets.wf.tf.targets $probtrack(output)/waypoints.txt
                    set flags "$flags --waypoints=$probtrack(output)/waypoints.txt "
	        } 
	        if { $probtrack(classify_yn) == 1 } {
                    fdt_exp w $w.data.targets.cf.tf.targets $probtrack(output)/targets.txt
                    set flags "$flags --targetmasks=$probtrack(output)/targets.txt --os2t "
                }

       		fdt_monitor_short $w "$FSLDIR/bin/probtrackx $flags"
                if { $probtrack(classify_yn) == 1 } {
	           fdt_monitor_short $w "$FSLDIR/bin/find_the_biggest ${logdir}/seeds_to_* biggest >> ${logdir}/fdt_seed_classification.txt"
		}
       	    }
	    if { $probtrack(mode) == "simple" } {
	        puts "rm ${filebase}_coordinates.txt"
	        exec rm ${filebase}_coordinates.txt
	    }
	    close $log
	    if { $copylog != "" } {
		puts "mv $logfile $copylog"
		exec mv $logfile $copylog
	    } else {
		puts "rm $logfile"
		exec rm $logfile
	    }
	}
	registration {
	    global registration

	    set errorStr ""
	    if { $registration(directory) == ""  } { set errorStr "You must specify the bedpost directory!" }
	    if { $registration(struct_yn) && $registration(struct_image) == ""  } { set errorStr "$errorStr You must specify the structural image!" }
	    if { $registration(standard_yn) && $registration(standard_image) == ""  } { set errorStr "$errorStr You must specify the standard image!" }
	    if { $errorStr != "" } {
		MxPause $errorStr
		return
	    }

	    exec mkdir -p [ file join $registration(directory) xfms ]
	    set eyefd [ open [ file join $registration(directory) xfms eye.mat ] w ]
	    puts $eyefd "1 0 0 0"
	    puts $eyefd "0 1 0 0"
	    puts $eyefd "0 0 1 0"
	    puts $eyefd "0 0 0 1"
	    close $eyefd

	    set diff2str   [ file join $registration(directory) xfms diff2str.mat ]
	    set str2diff   [ file join $registration(directory) xfms str2diff.mat ]
	    set str2stand  [ file join $registration(directory) xfms str2standard.mat ]
	    set stand2str  [ file join $registration(directory) xfms standard2str.mat ]
	    set diff2stand [ file join $registration(directory) xfms diff2standard.mat ]
	    set stand2diff [ file join $registration(directory) xfms standard2diff.mat ]
	    set diff       [ file join $registration(directory) nodif_brain ]
	    if { $registration(struct_yn) } {
		set searchrx  "-searchrx -$registration(struct_search) $registration(struct_search)"
		set searchry  "-searchry -$registration(struct_search) $registration(struct_search)"
		set searchrz  "-searchrz -$registration(struct_search) $registration(struct_search)"
		set options   "$searchrx $searchry $searchrz -dof $registration(struct_dof)"
		fdt_monitor $w "${FSLDIR}/bin/flirt -in $diff -ref $registration(struct_image) -omat $diff2str $options -cost $registration(struct_costfn)"
		fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $str2diff -inverse $diff2str"
		if { $registration(standard_yn) } {
		    set searchrx  "-searchrx -$registration(standard_search) $registration(standard_search)"
		    set searchry  "-searchry -$registration(standard_search) $registration(standard_search)"
		    set searchrz  "-searchrz -$registration(standard_search) $registration(standard_search)"
		    set options   "$searchrx $searchry $searchrz -dof $registration(standard_dof)"
		    fdt_monitor $w "${FSLDIR}/bin/flirt -in $registration(struct_image) -ref $registration(standard_image) -omat $str2stand $options -cost $registration(standard_costfn)"
		    fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $stand2str -inverse $str2stand"
		    fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $diff2stand -concat $str2stand $diff2str"
		    fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $stand2diff -inverse $diff2stand"
		}
	    } elseif { $registration(standard_yn) } {
		set searchrx  "-searchrx -$registration(standard_search) $registration(standard_search)"
		set searchry  "-searchry -$registration(standard_search) $registration(standard_search)"
		set searchrz  "-searchrz -$registration(standard_search) $registration(standard_search)"
		set options   "$searchrx $searchry $searchrz -dof $registration(standard_dof)"
		fdt_monitor $w "${FSLDIR}/bin/flirt -in $diff -ref $registration(standard_image) -omat $diff2stand $options"
		fdt_monitor $w "${FSLDIR}/bin/convert_xfm -omat $stand2diff -inverse $diff2stand"
	    }
	    puts "Done!"
	    # Fudge to make the logic work
	    set canwrite 1
	}
    }

    if { $canwrite } { 
	MxPause "  Done!  "
	update idletasks
    }

    if {$dialog == "destroy"} {
        fdt:destroy $w
    }
}


proc fdt:destroy { w } {
    destroy $w
}    

set debugging 0

while {[llength $argv] > 0 } {
    set flag [lindex $argv 0]
    switch -- $flag {
	"-debugging" {
	    set debugging 1
	    set argv [lrange $argv 1 end]
	    puts "Debug mode!"
	}
	default { break }
    }
}


wm withdraw .
if { [ info exists env(MRDATADIR) ] } {
    set MRDATADIR $env(MRDATADIR)
} else {
    set MRDATADIR ~/MRdata
}

fdt:dialog .fdt $argv
tkwait window .fdt
