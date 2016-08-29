# dimensionless box

reg machNumber          $MACHNUMBER

# physical constants
reg boltzmannConstant   1

# defined quantities
reg massComp            1 
reg domainSize          1 
reg ambientDens         1 
reg ambientTemp         1
reg adiabaticIndex      [expr {5.0/3.0}]        ;# mono-ionic gas (protons)
reg meanMolWeight       1                       ;# pure Hydrogen ion

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

reg maxStirTime         [expr { 10 * $turnTime}]
reg maxSimTime          [expr { 25 * $turnTime}]

reg dtinit              [expr {$turnTime / 1e6}]
reg dtmin               [expr {$turnTime / 1e10}]
reg dtmax               [expr {$turnTime / 1e2}] 

reg chkPtInterval       [expr {$maxSimTime / 400}]

# enforce scientifc notation with specified precision
proc scinot {quantity {prec 8}} {
    format "%-3.${prec}f" $quantity
}

proc perTurn {quantity} {
    global turnTime
    format {%5.2f t/t_c} [expr {$quantity/$turnTime}]
}

set report [subst {
> Physical Constants:
    Boltzmann Constant          [scinot $boltzmannConstant]

> Box Parameters
    domain size:                [scinot $domainSize]

> Hydrodynamics
    adiabatic index:            [format {%.3f} $adiabaticIndex]
    ambient temperature:        [scinot $ambientTemp]
    mass components:            [scinot $massComp]
    ambient density:            [scinot $ambientDens]
    ambient pressure:           [scinot $ambientPres]
    polytropic constant:        [scinot $polytropeConst]
    speed of sound:             [scinot $soundSpeed]

> Magnetodynamics:
    external magnetic field:
        Bx                      [scinot $Bx]
        By                      [scinot $By]
        Bz                      [scinot $Bz]

> Turbulent Stirring
    mach number                 [scinot $machNumber]
    characteristic velocity:    [scinot $maxVelocity]
    base forcing:               [scinot $baseForcing]
    turnover time:              [scinot $turnTime]

> Simulation
    chk point interval:         $chkPtInterval = [perTurn $chkPtInterval]
    max stirring time:          $maxStirTime = [perTurn $maxStirTime]
    max simulation time:        $maxSimTime = [perTurn $maxSimTime]
    
> Misc
    check point count:          [expr round($maxSimTime / $chkPtInterval)] 
}] 

# enforce scientifc notation with specified precision
proc scinot {quantity {prec 10}} {
    format "%.${prec}e" $quantity
}
