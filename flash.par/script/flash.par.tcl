#!/usr/bin/env tclsh

proc usage {} {
    puts "usage: [info script] <template file> <mach number> <definitions file> (json|machine)"
    exit 1
}

set NOW [clock seconds]

source utils.tcl
package require json

set TEMPLATEFILE    [utils::shift argv]
set MACHNUMBER      [utils::shift argv]
set DEFINITIONS     [utils::shift argv]
set OUTPUT          [utils::shift argv]

set meta [dict create]
proc reg {vname vdata} {
    global meta
    uplevel "set $vname [list $vdata]"
    dict set meta $vname $vdata
}

dict set meta {creation time} [clock format $NOW -format {%Y-%m-%d-%H-%M-%S}]
dict set meta cstamp $NOW

source $DEFINITIONS

switch -exact $OUTPUT {
    json    {puts [json::dict2json $meta]}
    machine { 
        puts {# _BEGIN_}
        dict for {k v} $meta {
            puts "# $k = $v"
        }
        puts {# _END_}
        puts {}
        puts "## [string repeat = 72] #"
        puts {}

        puts -nonewline [subst [utils::slurp $TEMPLATEFILE]]
    }
    default usage
}
