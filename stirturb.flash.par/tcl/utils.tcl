namespace eval utils {

# =========================================================================== #
# FILE AND DIRECTORY OPERATIONS #

proc mktmpdir {path} {
    set path [file normalize $path]
    append path /
    set chars abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789
    for {set i 0} {$i < 10} {incr i} {
        for {set j 0} {$j < 10} {incr j} {
            append path [string index $chars [expr {int(rand() * 62)}]]
        }
        if {![file exists $path]} {
            file mkdir $path
            file attributes $path -permissions 0700
            return $path
        }
    }
    error "failed to find an unused temporary directory name"
}

proc slurp {path} {
    set FH [open $path "r"]
    set data [read $FH]
    close $FH
    return $data
}

proc spurt {path data {mode w}} {
    set FH [open $path $mode]
    puts -nonewline $FH $data
    close $FH
    return 1
}

# =========================================================================== #
# LIST OPERATIONS #

proc map {body list} {
    lmap _ $list $body
}

proc grep {body list} {
    lmap _ $list {if [expr {$body}] {set _} {continue}}
}

proc head {list} {
    lindex $list 0
}

proc tail {list} {
    lreplace $list 0 0
}

# =========================================================================== #
# STACK OPERATIONS #

# remove first element
proc shift {lname} {
    upvar $lname mylist
    set fstel [lindex $mylist 0]
    set mylist [lreplace $mylist 0 0]
    return $fstel
}

# remove last element
proc pop {lname} {
    upvar $lname mylist
    set lastel [lindex $mylist end]
    set mylist [lreplace $mylist end end]
    return $lastel
}

# append element
proc push {lname el} {
    upvar $lname mylist
    set mylist [linsert $mylist end $el]
    return [llength $mylist]
}

# prepend element
proc pull {lname el} {
    upvar $lname mylist
    set mylist [linsert $mylist 0 $el]
    return [llength $mylist]
}

# cycle through list and return first element
proc cycle {lname {direction 1}} {
    upvar $lname mylist

    set el [head $mylist]
    # clock wise
    if {$direction == 1} {
        push mylist [shift mylist]
    # anti clock wise
    } elseif {$direction == -1} {
        pull mylist [pop mylist]
    }

    return $el
}

# =========================================================================== #
# LOGGER MECHANISMS #

proc dmsg {msg} {
    variable DEBUG_LEVEL
    if {[info exists DEBUG_LEVEL] && $DEBUG_LEVEL > 0} {
        puts stderr "\[DEBUG\] $msg"
    }
}

proc emsg {msg err} {
    puts stderr "\[ERROR\] \[[dict get $err -errorcode]\] $msg"
}


# =========================================================================== #
# CHARACTER ENCODING / DECODING OPERATIONS

proc utf8 {hex} {
    set hex [string map {% {}} $hex]
    encoding convertfrom utf-8 [binary decode hex $hex]
}

namespace eval url {
    variable map
    variable alphanumeric a-zA-Z0-9._~-
    namespace export encode decode
    namespace ensemble create
}

proc url::init {} {
        variable map
        variable alphanumeric a-zA-Z0-9._~-

        for {set i 0} {$i <= 256} {incr i} {
                set c [format %c $i]
                if {![string match \[$alphanumeric\] $c]} {
                        set map($c) %[format %.2x $i]
                }
        }
        # These are handled specially
        array set map { " " + \n %0d%0a }
}

url::init

proc url::encode {str} {
        variable map
        variable alphanumeric

        # The spec says: "non-alphanumeric characters are replaced by '%HH'"
        # 1 leave alphanumerics characters alone
        # 2 Convert every other character to an array lookup
        # 3 Escape constructs that are "special" to the tcl parser
        # 4 "subst" the result, doing all the array substitutions

        regsub -all \[^$alphanumeric\] $str {$map(&)} str
        # This quotes cases like $map([) or $map($) => $map(\[) ...
        regsub -all {[][{})\\]\)} $str {\\&} str
        return [subst -nocommand $str]
}

proc url::decode {str} {
        # rewrite "+" back to space
        # protect \ from quoting another '\'
        set str [string map [list + { } "\\" "\\\\"] $str]

        # prepare to process all %-escapes
        regsub -all -- {%([A-Fa-f0-9][A-Fa-f0-9])} $str {\\u00\1} str

        # process \u unicode mapped chars
        return [encoding convertfrom utf-8 [subst -novar -nocommand $str]]
}

proc main? {} {
    global argv0
    if {!([info exists argv0] && [file tail [info script]] eq [file tail $argv0])} {
        return 0
    } else {
        return 1
    }
}

proc set-srcpath {} {
    global SRCPATH
    set SRCPATH [file dirname [file normalize [info script]]]
}

proc rsplit {rx in} {
    set FS \x1c
    regsub -all $rx $in $FS out
    split $out $FS
}

}
