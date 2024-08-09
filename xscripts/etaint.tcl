# FG 20240807 v1.0

proc etaint { } {
    puts "ETA_INT> [tcloutr modkeyval eta_int]"
}

proc etainterr { args } {

    set niter 100
    if {[llength $args] > 0} {
       set niter [lindex $args 0]
    }
    set confidencelevel 68
    if {[llength $args] > 2} {
       set confidencelevel [lindex $args 2]
    }

    set savechatlevel [scan [tcloutr chatter] "%d"]
#    chatter 0

    # Get number of models and parameters
    set numpars [list]
    set mnames [tcloutr modnames]
    foreach mname $mnames {
        lappend numpars [tcloutr modpar $mname]
    }
    # Save current parameters to recover best-fit
    set i0 0
    set j 0
    foreach mname $mnames {
        for {set ipar 1} {$ipar <= [lindex $numpars $j]} {incr ipar} {
            scan [tcloutr param $mname:$ipar] "%f" sparval([expr {$i0+$ipar}])
            set sparislink([expr {$i0+$ipar}]) [string range [tcloutr plink $mname:$ipar] 0 0]
            set sparlink([expr {$i0+$ipar}]) [string trimleft [tcloutr plink $mname:$ipar] ?TF?]
        }
        incr j
        set i0 [expr {$i0+$ipar}]        
    }
  
#    set i0 0
#    set j 0   
#    foreach mname $mnames {
#        for {set ipar 1} {$ipar <= [lindex $numpars $j]} {incr ipar} {
#            puts $sparval([expr {$i0+$ipar}])
#            puts $sparislink([expr {$i0+$ipar}])
#            puts $sparlink([expr {$i0+$ipar}])
#        }
#        incr j
#        set i0 [expr {$i0+$ipar}]        
#    }    

    set etaints [list]
    set etaintavg 0.0

    for {set iter 1} {$iter <= $niter} { incr iter} {      
       chatter 0       
       tclout simpars

       regsub -all { +} [string trim $xspec_tclout] { } cpars
       set lpars [split $cpars]

       set i0 0
       set j 0
       foreach mname $mnames {
           set outstr " "
           for {set ipar 0} {$ipar < [lindex $numpars $j]} {incr ipar} {
               append outstr "& [lindex $lpars [expr {$ipar+$i0}]] "
           }
           newpar $mname:1-[lindex $numpars $j] $outstr
           incr j
           set i0 [expr {$i0+$ipar}]
       }

       set etaint [tcloutr modkeyval eta_int]       
       set etaintavg [expr $etaintavg+$etaint/$niter]             
       lappend etaints $etaint
       
       puts "ETA_INT> $iter of $niter: $etaint"

       set i0 0
       set j 0
       foreach mname $mnames {
           set outstr " "
           for {set ipar 1} {$ipar <= [lindex $numpars $j]} {incr ipar} {
               if {$sparislink([expr {$ipar+$i0}]) == "F"} {
                   append outstr "& $sparval([expr {$ipar+$i0}]) "
               } else {
                   append outstr "& $sparlink([expr {$ipar+$i0}]) "
               }
           }
           newpar $mname:1-[lindex $numpars $j] $outstr
           incr j
           set i0 [expr {$i0+$ipar}]
       }     
       fit
       chatter $savechatlevel
    }
    
    set sorted_etaints [lsort -real -increasing $etaints]        
    
    set low [expr ($niter*(100-$confidencelevel)/200)-1]
    set high [expr $niter-$low-1]
    set med [expr int(0.5*$niter)]

    puts " "
    puts "Median ($confidencelevel% C.I.) on ETA_INT : "
    puts "   [lindex $sorted_etaints $med] ([lindex $sorted_etaints $low], [lindex $sorted_etaints $high])"
    puts " <ETA_INT> = $etaintavg"
    puts " "


#    chatter $savechatlevel
   
}
