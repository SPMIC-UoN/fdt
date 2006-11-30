#   FSL interface for FDT (BEDPOST and ProbTrack)
#
#   Timothy Behrens, Heidi Johansen-Berg, Dave Flitney and Matthew Webster FMRIB Image Analysis Group
#
#   Copyright (C) 2006 University of Oxford
#
#   TCLCOPYRIGHT

#TO DO replace -filetypes * with -filetypes { } for directory selectors
source [ file dirname [ info script ] ]/fslstart.tcl
set TCLPATH [file dirname [ info script ] ]
regsub tcl $TCLPATH bin BINPATH
regsub tcl $TCLPATH doc/fdt HTMLPATH

set VERSION "1.0a"

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

    set LWIDTH 27
    set PWD [ pwd ]
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
    optionMenu2 $w.tool.menu probtrack(tool) -command "fdt:select_tool $w" eddy_current "Eddy current correction" bedpost "BEDPOST Estimation of diffusion parameters"  registration "Registration" probtrack "ProbTrack Probabilistic tracking" xutilssx "----------------------------------------------------" dtifit "DTIFit Reconstruct diffusion tensors" 
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

    FileEntry $w.registration.directory -textvariable registration(directory) -label "BEDPOST directory:" -title "Choose directory"  -width 35 -filetypes * -command "registration_set_directory $w" 

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

    FileEntry $w.ecc.input -textvariable eddy(input) -label "Diffusion weighted data:" -title "Choose diffusion weighted image"  -width 35 -filetypes IMAGE -command "ecc_update_files $w"

    FileEntry $w.ecc.output -textvariable eddy(output) 	-label "Corrected output data:" -title  "Choose output image name" -width 35 -filetypes IMAGE -command "ecc_update_files $w"

   set eddy(refnum) 0
   LabelSpinBox  $w.ecc.refnum -label "Reference volume"  -textvariable eddy(refnum) -range {0 100 1 } -width 6 

    pack $w.ecc.input $w.ecc.output $w.ecc.refnum -side top -padx 3 -pady 3 -expand yes -anchor w

    #------- DTIFit --------

    frame $w.dtifit

    FileEntry $w.dtifit.directory -textvariable dtifit(directory) -label  "Input directory:" -title "Choose directory" -width 35 -command "set_working_directory dtifit(cwd)"

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
	upvar $cwd myCWD

	set dirname [file dirname $filename]
	puts "switching from $myCWD to $dirname" 
	set myCWD $dirname
    }

    proc dtifit_update_files { w filename } {
	global dtifit

	set dtifit(output) [ file join [file dirname $dtifit(input)] dti ]
	set_working_directory dtifit(cwd) $dtifit(input)
    }
    
    set dtifit(cwd) $PWD

#All the below orignally had -directory $dtifit(cwd) 
    FileEntry $w.dtifit.expert.input -textvariable dtifit(input) -label  "Diffusion weighted data:" -title "Choose diffusion weighted image" -width 35 -filetypes IMAGE -command "dtifit_update_files $w" 
    $w.dtifit.expert.input.labf configure -width $LWIDTH

    FileEntry $w.dtifit.expert.mask -textvariable dtifit(mask) -label "BET binary brain mask:" -title "Choose BET brain mask file" -width 35 -filetypes IMAGE -command "set_working_directory dtifit(cwd)"
    $w.dtifit.expert.mask.labf configure -width $LWIDTH
    FileEntry $w.dtifit.expert.output -textvariable dtifit(output) -label "Output basename:" -title  "Choose output base name" -width 35 -filetypes * -command "set_working_directory dtifit(cwd)"
    $w.dtifit.expert.output.labf configure -width $LWIDTH
    FileEntry $w.dtifit.expert.bvecs -textvariable dtifit(bvecs) -label "Gradient directions:" -title  "Choose bvecs file" -width 35 -filetypes * -command "set_working_directory dtifit(cwd)"
    $w.dtifit.expert.bvecs.labf configure -width $LWIDTH
    FileEntry $w.dtifit.expert.bvals -textvariable dtifit(bvals) -label  "b values:" -title  "Choose bvals file" -width 35 -filetypes * -command "set_working_directory dtifit(cwd)"
    $w.dtifit.expert.bvals.labf configure -width $LWIDTH        

    pack $w.dtifit.expert.input $w.dtifit.expert.mask $w.dtifit.expert.output \
	$w.dtifit.expert.bvecs $w.dtifit.expert.bvals \
	-in $w.dtifit.expert -side top -padx 3 -pady 3 -expand yes -anchor w

    pack $w.dtifit.directory $w.dtifit.expert_yn -in $w.dtifit \
	-side top -padx 3 -pady 3 -expand yes -anchor w

    #------- BEDPOST --------

    frame $w.bedpost

    FileEntry $w.bedpost.directory -textvariable bedpost(directory) -label "Input directory:" -title "Choose directory" -width 35 -filetypes * -command "set_working_directory dtifit(cwd)"

    set bedpost(ecc_yn) 0

    pack $w.bedpost.directory \
	-in $w.bedpost -side top -padx 3 -pady 3 -expand yes -anchor w

    #-------- ...ProbTrack... -------- 

    NoteBook $w.probtrack -bd 2 -tabpady {5 10} -arcradius 3
    $w.probtrack insert 0 data -text "Data"
    $w.probtrack insert 1 options -text "Options"

    #-------- ...Mode specific options... --------

    set dataf [$w.probtrack getframe data]
    frame $w.data

    LabelFrame  $w.data.mode -text "Mode: "

    if { [ file exists [ file join $FSLDIR bin reord_OM ] ] } {
        optionMenu2 $w.data.mode.menu probtrack(mode) -command "fdt:probtrack_mode $w" xtitlex "Path distribution estimation" simple "  Single seed voxel" all "  Seed mask" maska "  Seed mask and waypoint masks" masks "  Two masks - symmetric" xutilssx "--------------------------------------------" seeds "Connectivity-based seed classification" xmatrixx "--------------------------------------------" mat1 "Matrix 1" mat2 "Matrix 2" mskmat "Mask matrix"
      $w.data.mode.menu.menu entryconfigure 0 -state disabled -background black
      $w.data.mode.menu.menu entryconfigure 5 -state disabled -background black
      $w.data.mode.menu.menu entryconfigure 7 -state disabled -background black
    } else {
    optionMenu2 $w.data.mode.menu probtrack(mode) -command "fdt:probtrack_mode $w" xtitlex "Path distribution estimation" simple "  Single seed voxel" all "  Seed mask" maska "  Seed mask and waypoint masks" masks "  Two masks - symmetric" xutilssx  "--------------------------------------------" seeds "Connectivity-based seed classification"
     
     $w.data.mode.menu.menu entryconfigure 0 -state disabled -background black
     $w.data.mode.menu.menu entryconfigure 5 -state disabled -background black

    }
    pack $w.data.mode.menu

    set probtrack(xfm) ""
    set probtrack(basename) "merged"
    set probtrack(mask) "nodif_brain_mask"

    proc probtrack_update_files { w filename } {
	global probtrack
	global FSLDIR

	if { ($probtrack(bedpost_dir) != "") && ($probtrack(seed) != "") } {
	    set probtrack(dir) \
		[ file join $probtrack(bedpost_dir) [ file tail [ exec $FSLDIR/bin/remove_ext $probtrack(seed) ] ] ]
	}
    }

    FileEntry $w.data.directory -textvariable probtrack(bedpost_dir) -label "BEDPOST directory" -title "Choose BEDPOST directory" -width 35 -filetypes * -command "probtrack_update_files $w"

    TitleFrame $w.data.target -text "Target list"    
    listbox $w.data.targets -height 6 -width 40 -yscrollcommand "$w.data.sb set"
    scrollbar $w.data.sb -command "$w.data.targets yview " 
    frame $w.data.tb
    button $w.data.tb.add -text "Add..."  -command "feat_file:setup_dialog $w a a a [namespace current] IMAGE {Select File} {fdt_add $w $w.data.targets} {}" 
    button $w.data.tb.del -text "Remove"  -command "fdt_sub $w $w.data.targets" 
    button $w.data.tb.imp -text "Load list..." -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_imp $w $w.data.targets} {}"
    button $w.data.tb.exp -text "Save list..." -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_exp $w $w.data.targets} {}"
    pack $w.data.tb.add $w.data.tb.del $w.data.tb.imp $w.data.tb.exp -side left
    pack $w.data.tb -in [$w.data.target getframe ] -side bottom  -expand yes -fill x -anchor w -padx 3 -pady 3
    pack $w.data.targets $w.data.sb -in [$w.data.target getframe ] -side left  -expand yes -fill y -anchor w -padx 3 -pady 3

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

    TitleFrame $w.data.seedspace -text "Seed space"
    FileEntry $w.data.seed -textvariable probtrack(seed) -label "Seed image:" -title "Choose seed image" -width 35 -filetypes IMAGE -command "probtrack_update_files $w"
    FileEntry $w.data.seed2 -textvariable probtrack(seed2) -label "Target image:" -title "Choose target image" -width 35 -filetypes IMAGE -command "probtrack_update_files $w"
 
    TitleFrame $w.data.waypoint -text "Waypoint masks"
    listbox $w.data.waypoints -height 6 -width 40 -yscrollcommand "$w.data.waypoint.sb set"
    scrollbar $w.data.waypoint.sb -command "$w.data.waypoints yview "
    frame $w.data.waybut
    button $w.data.waybut.add -text "Add..."  -command "feat_file:setup_dialog $w a a a [namespace current] IMAGE {Select File} {fdt_add $w $w.data.waypoints} {}" 
    button $w.data.waybut.del -text "Remove"  -command "fdt_sub $w $w.data.waypoints" 
    button $w.data.waybut.imp -text "Load list..." -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_imp $w $w.data.waypoints} {}"
    button $w.data.waybut.exp -text "Save list..." -command "feat_file:setup_dialog $w a a a [namespace current] * {Select File} {fdt_exp $w $w.data.waypoints} {}" 
    pack $w.data.waybut.add $w.data.waybut.del $w.data.waybut.imp $w.data.waybut.exp -in $w.data.waybut -side left
    pack $w.data.waybut -in [$w.data.waypoint getframe ] -side bottom -expand yes -fill x -anchor w -padx 3 -pady 3
    pack $w.data.waypoints $w.data.waypoint.sb -in [$w.data.waypoint getframe ] -side left -expand yes -fill y -anchor w -padx 3 -pady 3

    set probtrack(x) 0
    set probtrack(y) 0
    set probtrack(z) 0
    set probtrack(units) vox

    #Co-ordinate edit frame
    LabelFrame   $w.data.seedxyz -text "Seed:"
    LabelSpinBox $w.data.seedxyz.x -label "X" -textvariable probtrack(x) -range {-1000000 1000000 1 } -width 6
    LabelSpinBox $w.data.seedxyz.y -label "Y" -textvariable probtrack(x) -range {-1000000 1000000 1 } -width 6
    LabelSpinBox $w.data.seedxyz.z -label "Z" -textvariable probtrack(x) -range {-1000000 1000000 1 } -width 6
    radiobutton $w.data.seedxyz.vox -text "vox" -value vox -variable probtrack(units)
    radiobutton $w.data.seedxyz.mm  -text "mm"  -value mm  -variable probtrack(units)
    pack $w.data.seedxyz.x $w.data.seedxyz.y $w.data.seedxyz.z $w.data.seedxyz.vox $w.data.seedxyz.mm -side left -padx 2 -pady 0

    proc probtrack_toggle_reference { w } {
	global probtrack

	if { $probtrack(usereference_yn) } { 
	    if { $probtrack(mode) == "simple" } {
		pack $w.data.reference $w.data.xfm -in [$w.data.seedspace getframe ] \
		    -side top -after $w.data.usereference_yn -padx 3 -pady 3
	    } else {
		pack $w.data.xfm -in [$w.data.seedspace getframe ] \
		    -side top -after $w.data.usereference_yn -padx 3 -pady 3
	    }
	} else { 
	    pack forget $w.data.reference
	    pack forget $w.data.xfm
	}
    }

    checkbutton $w.data.usereference_yn -text "Seed space is not diffusion" \
	-variable probtrack(usereference_yn) \
	-command "probtrack_toggle_reference $w"

    FileEntry $w.data.reference -textvariable probtrack(reference) -label "Seed space reference image:" -title "Choose reference image" -width 35 -filetypes IMAGE 

    FileEntry $w.data.xfm -textvariable probtrack(xfm) -label "Seed to diff-space transform:" -title "Select seed-space to DTI-space transformation matrix" -width 35 -filetypes *.mat

    proc probtrack_toggle_exclude { w } {
	global probtrack

	if { $probtrack(exclude_yn) } { 
	    pack $w.data.exclude -in [$w.data.seedspace getframe ] \
		-side top -after $w.data.exclude_yn -padx 3 -pady 3
	} else { 
	    pack forget $w.data.exclude
	}
    }

    checkbutton $w.data.exclude_yn -text "Use exclusion mask" -variable probtrack(exclude_yn) \
	-command "probtrack_toggle_exclude $w"

    FileEntry $w.data.exclude -textvariable probtrack(exclude) -label "Exclusion mask:" -title "Select exclusion mask image" -width 35 -filetypes IMAGE 

    pack $w.data.seed $w.data.usereference_yn $w.data.exclude_yn -in\
	[$w.data.seedspace getframe ] -padx 3 -pady 3 -side top -anchor w

    probtrack_toggle_reference $w
    probtrack_toggle_exclude $w

    FileEntry $w.data.lrmask -textvariable probtrack(lrmask) -label "Low resolution mask:" -title  "Choose low resolution mask" -width 35 -filetypes IMAGE 

    pack $w.data.lrmask \
	-in [$w.data.seedspace getframe ] \
	-side top -padx 3 -pady 3 -anchor w

    TitleFrame $w.data.output -text "Ouput"

   FileEntry $w.data.dir -textvariable probtrack(dir) -label  "Output directory:" -title  "Name the output directory" -width 35 -filetypes *

   FileEntry $w.data.out -textvariable probtrack(output) -label "Output:" -title "Choose output file name" -width 35 -filetypes IMAGE

    pack $w.data.out $w.data.dir -in [$w.data.output getframe ] -side top -expand yes -fill x -padx 3 -pady 3 -anchor w

    pack $w.data.mode $w.data.directory $w.data.target $w.data.seedspace $w.data.output \
	-in $w.data -side top -padx 3 -pady 3 -anchor nw

    pack $w.data -in $dataf -side top -padx 3 -pady 3 -anchor nw -expand yes -fill both

    #-------- ...Options... --------

    set optionsf [$w.probtrack getframe options]

    TitleFrame $w.options -text "Basic Options"

    checkbutton $w.options.verbose \
	    -text "Verbose" \
	    -variable probtrack(verbose_yn)
    
    set probtrack(nparticles) 5000
    LabelSpinBox $w.options.nparticles -label  "Number of samples" -textvariable probtrack(nparticles) -range {1 1e24 100 } -width 6 

    set probtrack(curvature) 0.2
    LabelSpinBox $w.options.curvature -label "Curvature threshold" -textvariable probtrack(curvature) -range {0.0 1.0 0.01 } -width 4

    set probtrack(loopcheck_yn) 1
    checkbutton $w.options.loopcheck \
	-text "Loopcheck" \
	-variable probtrack(loopcheck_yn)

    collapsible frame $w.advanced -title "Advanced Options"

    set probtrack(nsteps) 2000
    LabelSpinBox $w.advanced.nsteps -label "Maximum number of steps" -textvariable probtrack(nsteps) -range {2 1000000 10 } -width 6

    set probtrack(steplength) 0.5
    LabelSpinBox $w.advanced.steplength -label "Step length (mm)" -textvariable probtrack(steplength) -range {0 10000 0.1} -width 4

    set probtrack(modeuler_yn) 0
    checkbutton $w.advanced.modeuler \
	-text "Use modified Euler streamlining" \
	-variable probtrack(modeuler_yn)

    pack \
	$w.advanced.modeuler \
	$w.advanced.nsteps \
	$w.advanced.steplength \
	-in $w.advanced.b -side top -pady 3 -padx 6 -anchor nw

    pack \
	$w.options.nparticles \
	$w.options.curvature \
	$w.options.verbose \
	$w.options.loopcheck \
	-in [$w.options getframe ] -side top -pady 3 -padx 6 -anchor nw

    pack $w.options $w.advanced -in $optionsf -side top -pady 3 -padx 6 -anchor nw -expand yes -fill both

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

    set seedspacef [$w.data.seedspace getframe ]

    pack forget $w.data.target
    pack forget $w.data.lrmask
    pack forget $w.data.usereference_yn
    pack forget $w.data.reference
    pack forget $w.data.seed2
    pack forget $w.data.waypoint
    pack forget $w.data.seedxyz
    pack forget $w.data.out
    pack $w.data.seed $w.data.usereference_yn $w.data.xfm -in $seedspacef -side top \
	-before $w.data.exclude_yn -padx 3 -pady 3 -anchor nw
    $w.data.seed configure -label "Seed image:"
    pack $w.data.dir -in [$w.data.output getframe ] -side top -padx 3 -pady 3 -anchor nw
    switch -- $probtrack(mode) {
  	simple {
 	    pack $w.data.seedxyz -in $seedspacef -side top -padx 3 -pady 3 -before $w.data.usereference_yn
	    pack $w.data.reference -in $seedspacef -side top -padx 3 -pady 3 -before $w.data.xfm -anchor nw
	    pack $w.data.out -in [$w.data.output getframe] -side top -padx 3 -pady 3 -anchor nw
  	    pack forget $w.data.dir
	    pack forget $w.data.seed
  	}
	maska {
	    pack $w.data.waypoint -in $w.data -side top -padx 3 -pady 3 -after $w.data.seedspace -anchor nw
	}
	masks {
	    pack $w.data.seed2 -in $seedspacef -side top -after $w.data.seed
	    $w.data.seed  configure -label "Mask image 1:"
	    $w.data.seed2 configure -label "Mask image 2:"
	}
  	seeds {
  	    pack $w.data.target -in $w.data -side top -padx 3 -pady 3 -after $w.data.seedspace -anchor nw
  	}
  	mat2 {
  	    pack $w.data.lrmask -in $seedspacef \
 		-side top -padx 3 -pady 3 -after $w.data.xfm -anchor nw
  	}
    }
    probtrack_toggle_reference $w
    probtrack_toggle_exclude $w
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

    global probtrack(tool)
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
		puts "bedpost $bedpost(directory)"
		fdt_monitor $w "${FSLDIR}/bin/bedpost $bedpost(directory)"
	    }
	}
	probtrack {
	    global probtrack

	    set errorStr ""
	    if { $probtrack(bedpost_dir) == ""  } { set errorStr "You must specify the bedpost directory!" }
	    if { $probtrack(usereference_yn) && $probtrack(xfm) == "" } { set errorStr "$errorStr You must specify the seed-diff transform!" }
	    if { $probtrack(exclude_yn) && $probtrack(exclude) == "" } { set errorStr "$errorStr You must specify the exclusion mask!" }

	    set flags ""
	    if { $probtrack(verbose_yn) == 1 } { set flags "$flags -V 1" }
	    if { $probtrack(loopcheck_yn) == 1 } { set flags "$flags -l" }
#	    if { $probtrack(usef_yn) == 1 } { set flags "$flags -f" }
	    if { $probtrack(modeuler_yn) == 1 } { set flags "$flags --modeuler" }
	    set flags "$flags -c $probtrack(curvature) -S $probtrack(nsteps) --steplength=$probtrack(steplength) -P $probtrack(nparticles)"

	    set tn [open "| $BINPATH/tmpnam"]
	    gets $tn filebase
	    close $tn
	    set logfile "${filebase}_log.tcl"
	    set log [open "$logfile" w]
	    puts $log "set tool $probtrack(tool)"
	    set copylog ""
	    set ssopts ""
	    if { $probtrack(usereference_yn) } {
		set ssopts "--xfm=$probtrack(xfm)"
	    }
	    if { $probtrack(exclude_yn) == 1 } {
		set ssopts "$ssopts --rubbish=$probtrack(exclude)"
		puts $log "set probtrack(exclude) $probtrack(exclude)"
	    }
	    set basics "--forcedir -s $probtrack(bedpost_dir)/merged -m $probtrack(bedpost_dir)/nodif_brain_mask"	    

    	    foreach entry {bedpost_dir xfm mode exclude_yn usereference_yn verbose_yn loopcheck_yn modeuler_yn \
			       curvature nsteps steplength nparticles} {
		puts $log "set probtrack($entry) $probtrack($entry)"
	    }

	    switch -- $probtrack(mode) {
		simple {
		    if { $probtrack(output) == ""  } { set errorStr "$errorStr You must specify the output basename!" }
		    if { $probtrack(usereference_yn) && $probtrack(reference) == "" } { 
			set errorStr "$errorStr You must specify a seed space reference image!" 
		    }
		    if { $probtrack(usereference_yn) == 1 } {
			set ssopts "--xfm=$probtrack(xfm) --seedref=$probtrack(reference)"
			puts $log "set probtrack(reference) $probtrack(reference)"
			puts $log "set probtrack(xfm) $probtrack(xfm)"
		    } else {
			set ssopts ""
		    }
		    if { $probtrack(exclude_yn) == 1 } {
			set ssopts "$ssopts --rubbish=$probtrack(exclude)"
		    }
		    set fd [ open "${filebase}_coordinates.txt" w ]
		    $w.data.seedxyz update
		    if { $probtrack(units) == "mm" } {
			set x $probtrack(x)
			set y $probtrack(y)
			set z $probtrack(z)
			if { $probtrack(usereference_yn) && $probtrack(reference) != "" } {
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
		    puts $log "set probtrack(output) $probtrack(output)"

		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }

		    set canwrite 1
		    if { [ file exists $probtrack(output) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(output)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(output)"
			    exec rm -rf $probtrack(output)
			}
		    }
		    if { $canwrite } {
			set copylog "fdt.log"
			fdt_monitor_short $w "$FSLDIR/bin/probtrack --mode=simple -x ${filebase}_coordinates.txt $basics $ssopts $flags -o $probtrack(output)"
		    }
		    puts "rm ${filebase}_coordinates.txt"
		    exec rm ${filebase}_coordinates.txt
		}
		all {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=seedmask -x $probtrack(seed) $basics $ssopts $flags -o fdt_paths --dir=$probtrack(dir)"
		    }
		}
		masks {
		    if { $probtrack(seed2) == "" } { set errorStr "$errorStr You must select a target image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(seed2) $probtrack(seed2)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=twomasks_symm --seed=$probtrack(seed) --mask2=$probtrack(seed2) $basics $ssopts $flags -o fdt_paths --dir=$probtrack(dir)"
		    }
		}
		maska {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			puts $log "$w.data.watpoints load \"$probtrack(dir)/waypoints.txt\""
			$w.data.waypoints save "$probtrack(dir)/waypoints.txt"
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=waypoints --seed=$probtrack(seed) --mask2=${probtrack(dir)}/waypoints.txt $basics $ssopts $flags -o fdt_paths --dir=$probtrack(dir)"
		    }
		}
		seeds {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"


		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			puts $log "$w.data.targets load \"$probtrack(dir)/targets.txt\""
			$w.data.targets save "$probtrack(dir)/targets.txt"
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=seeds_to_targets -x $probtrack(seed) $basics $ssopts $flags --targetmasks=${probtrack(dir)}/targets.txt --dir=$probtrack(dir)"
		    }
		}
		mat1 {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=matrix1 -x $probtrack(seed) $basics $ssopts $flags -o fdt_matrix --dir=$probtrack(dir)"
		    }
		}
		mat2 {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $probtrack(lrmask) == "" } { set errorStr "$errorStr You must specify the low resolution mask!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"
		    puts $log "set probtrack(lrmask) $probtrack(lrmask)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=matrix2 -x $probtrack(seed) $basics $ssopts --lrmask=$probtrack(lrmask) $flags -o fdt_matrix --dir=$probtrack(dir)"
		    }
		}
		mskmat {
		    if { $probtrack(seed) == ""  } { set errorStr "$errorStr You must specify a seed image!" }
		    if { $probtrack(dir) == ""  } { set errorStr "$errorStr You must specify the output directory!" }
		    if { $errorStr != "" } {
			MxPause $errorStr
			return
		    }
		    puts $log "set probtrack(seed) $probtrack(seed)"
		    puts $log "set probtrack(dir) $probtrack(dir)"

		    set canwrite 1
		    if { [ file exists $probtrack(dir) ] } {
			set canwrite [  YesNoWidget "Overwrite $probtrack(dir)?" Yes No ]
			if { $canwrite } {
			    puts "rm -rf $probtrack(dir)"
			    exec rm -rf $probtrack(dir)
			}
		    }
		    if { $canwrite } {
			puts "mkdir -p $probtrack(dir)"
			exec mkdir -p $probtrack(dir)
			set copylog "$probtrack(dir)/fdt.log"
			fdt_monitor $w "$FSLDIR/bin/probtrack --mode=maskmatrix -x $probtrack(seed) $basics $ssopts $flags -o fdt_matrix --dir=$probtrack(dir)"
		    }
		}
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
