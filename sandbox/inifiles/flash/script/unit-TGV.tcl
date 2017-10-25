# dimensionless box

reg xmin  0.0
reg xmax  1.0
reg domsize [expr {$xmax - $xmin}]

reg mu0   1.0
reg mach  $MACHNUMBER
reg gamma 1.4

reg dens0 1.0
reg pres0 1.0 
reg temp0 1.0

reg sona0 [expr {sqrt($pres0/$dens0)}]

reg velmax [expr {$mach * $sona0}]

reg velx0  $velmax 
reg vely0 -$velmax
reg velz0 0.0

reg accx0 0.0
reg accy0 0.0
reg accz0 0.0

reg magx0 0.0
reg magy0 0.0
reg magz0 0.0

# derived quantities

reg dyntime             [expr {$domsize / $velmax}]

reg maxSimTime          [expr {6 * $dyntime}]

reg dtinit              [expr {$dyntime / 1e6}]
reg dtmin               [expr {$dyntime / 1e10}]
reg dtmax               [expr {$dyntime / 1e2}] 

reg chkPtInterval       [expr {$maxSimTime / 100}]

# forcing

reg maxStirTime         [expr {4 * $dyntime}]
reg baseforcing         [expr {$velmax**2 / $domsize}]

# enforce scientifc notation with specified precision
proc scinot {quantity {prec 10}} {
    format "%.${prec}e" $quantity
}
