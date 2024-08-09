# FG 20240807 v1.0

proc etaint12 { } {
    puts "ETA_INT12> (1) [tcloutr modkeyval eta_int1]  (2) [tcloutr modkeyval eta_int2]"
}

proc etaint12err { args } {

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

    set etaints1 [list]
    set etaints2 [list]        

    set etaint1avg 0.0
    set etaint2avg 0.0    

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
       
       set etaint1 [tcloutr modkeyval eta1_int]
       set etaint2 [tcloutr modkeyval eta2_int]       
       set etaint1avg [expr $etaint1avg+$etaint1/$niter]
       set etaint2avg [expr $etaint2avg+$etaint2/$niter]             
       lappend etaints1 $etaint1
       lappend etaints2 $etaint2           

       puts "ETA_INT12> $iter of $niter: (1) $etaint1  (2) $etaint2"

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
    
    set sorted_etaints1 [lsort -real -increasing $etaints1]
    set sorted_etaints2 [lsort -real -increasing $etaints2]        
    
    set low [expr ($niter*(100-$confidencelevel)/200)-1]
    set high [expr $niter-$low-1]
    set med [expr int(0.5*$niter)]

    puts " "
    puts "Median ($confidencelevel% C.I.) on ETA_INT12 : "
    puts "(1)   [lindex $sorted_etaints1 $med] ([lindex $sorted_etaints1 $low], [lindex $sorted_etaints1 $high])"
    puts " <ETA_INT1> = $etaint1avg"
    puts "(2)   [lindex $sorted_etaints2 $med] ([lindex $sorted_etaints2 $low], [lindex $sorted_etaints2 $high])"
    puts " <ETA_INT2> = $etaint2avg"    
    puts " "


#    chatter $savechatlevel
   
}
