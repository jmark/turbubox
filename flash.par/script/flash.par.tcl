#!/usr/bin/env tclsh

proc usage {} {
    puts "usage: [info script] <template file> <mach number> (unit|cgs) (human|json|machine)"
    exit 1
}

set NOW [clock seconds]

source utils.tcl
package require json

set TEMPLATEFILE    [utils::shift argv]
set MACHNUMBER      [utils::shift argv]
set UNITSYSTEM      [utils::shift argv]
set OUTPUT          [utils::shift argv]

set meta [dict create]
proc reg {vname vdata} {
    global meta
    uplevel "set $vname [list $vdata]"
    dict set meta $vname $vdata
}

## ... params + report

switch -exact $UNITSYSTEM {
    unit    {source unit.tm}
    cgs     {source cgs.tm}
    default usage 
}

switch -exact $OUTPUT {
    human   {puts $report}
    json    {puts [json::dict2json $meta]}
    machine { 
        puts "#ctime=[clock format $NOW -format {%Y-%m-%d-%H-%M-%S}]"
        puts "#cstamp=$NOW"
        puts "#unitsystem=$UNITSYSTEM"

        dict for {k v} $meta {
            puts "#$k=$v"
        }
        puts {}

        puts "## [string repeat = 72] #"
        puts {## Human Readable Report}
        puts "## [string repeat = 72] #"
        puts -nonewline {##}
        puts [join [split $report "\n"] "\n## "]
        puts "## [string repeat = 72] #"
        puts {}

        puts -nonewline [subst [utils::slurp $TEMPLATEFILE]]
    }
    default usage
}
