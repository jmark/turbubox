#!/usr/bin/env tclsh

source utils.tcl

package require json

set meta [dict create]
proc reg {vname vdata} {
    global meta
    uplevel "set $vname [list $vdata]"
    dict set meta $vname $vdata
}

## We use the CGS-Units System

# C ->  Centimeter   (Length)
# G ->  Gramm        (Mass)
# S ->  Second       (Time)
#       Kelvin       (Temperature)
#       Erg          (Energy: 10e-7 Joules)

reg machNumber          [lindex $argv 1]

# physical constants
reg boltzmannConstant   1.38064852e-16  ;# [erg/K]
reg massProton          1.672621898e-24 ;# [g]
reg parsec              3.0857e18       ;# [cm]

# defined quantities
reg massComp            $massProton             ;# mass of the gas components
reg domainSize          [expr {10 * $parsec}]   ;# cubic periodic box
reg ambientDens         [expr {10 * $massComp}] ;# 10 protons per ccm 
reg ambientTemp         15                      ;# Kelvin
reg adiabaticIndex      [expr {5.0/3.0}]        ;# mono-ionic gas (protons)
reg meanMolWeight       1                       ;# pure Hydrogen ion

# external magnetic field flux density
reg Bx 0.0
reg By 0.0
reg Bz 0.0

# derived quantities
reg polytropeConst      [expr {$boltzmannConstant * $ambientTemp / $massComp}] 
reg ambientPres         [expr {$polytropeConst * $ambientDens}] 
reg soundSpeed          [expr {sqrt($polytropeConst)}]

reg maxVelocity         [expr {$machNumber * $soundSpeed}]
reg turnTime            [expr {$domainSize / $maxVelocity}]

reg baseForcing         [expr {$maxVelocity**2 / $domainSize}]

reg maxStirTime         [expr {  2 * $turnTime}]
reg maxSimTime          [expr { 10 * $turnTime}]

reg dtinit              [expr {$turnTime / 1e6}]
reg dtmin               [expr {$turnTime / 1e10}]
reg dtmax               [expr {$turnTime / 1e2}] 

reg chkPtInterval       [expr {$maxSimTime / 400}]

# enforce scientifc notation with specified precision
proc scinot {quantity {prec 3}} {
    format "%-3.${prec}e" $quantity
}

proc Myr {seconds} {
    format {%6.2f Myr} [expr {$seconds / (60.0 * 60 * 24 * 365) / 1e6}]
}

proc perTurn {quantity} {
    global turnTime
    format {%5.2f t/t_c} [expr {$quantity/$turnTime}]
}

if {[lindex $argv 2] eq {-json}} {
    puts stderr [json::dict2json $meta]
} else {
puts stderr [subst {
> Physical Constants:
    Boltzmann Constant          [scinot $boltzmannConstant] erg/K

> Box Parameters
    domain size:                [scinot $domainSize] cm = [expr {$domainSize / $parsec}] pc

> Hydrodynamics
    adiabatic index:            [format {%.3f} $adiabaticIndex]
    ambient temperature:        $ambientTemp K
    mass components:            [scinot $massComp] g
    ambient density:            [scinot $ambientDens] g/ccm 
    ambient pressure:           [scinot $ambientPres] Ba
    polytropic constant:        [scinot $polytropeConst] erg/g
    speed of sound:             [scinot $soundSpeed] cm/s

> Magnetodynamics:
    external magnetic field:
        Bx                      [scinot $Bx] G
        By                      [scinot $By] G
        Bz                      [scinot $Bz] G

> Turbulent Stirring
    mach number                 $machNumber
    characteristic velocity:    [scinot $maxVelocity] cm/s
    base forcing:               [scinot $baseForcing] cm/s^2
    turnover time:              [scinot $turnTime] s = [Myr $turnTime]

> Simulation
    chk point interval:         [Myr $chkPtInterval] = [perTurn $chkPtInterval]
    max stirring time:          [Myr $maxStirTime] = [perTurn $maxStirTime]
    max simulation time:        [Myr $maxSimTime] = [perTurn $maxSimTime]
    
> Misc
    check point count:          [expr round($maxSimTime / $chkPtInterval)] 
}] 
}

# enforce scientifc notation with specified precision
proc scinot {quantity {prec 10}} {
    format "%.${prec}e" $quantity
}

puts -nonewline [subst [utils::slurp [lindex $argv 0]]]
