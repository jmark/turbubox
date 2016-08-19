import math

## We use the CGS-Units System

# C ->  Centimeter   (Length)
# G ->  Gramm        (Mass)
# S ->  Second       (Time)
#       Kelvin       (Temperature)
#       Erg          (Energy: 10e-7 Joules)

# physical constants
boltzmannConstant   = 1.38064852e-16  # [erg/K]
massProton          = 1.672621898e-24 # [g]
parsec              = 3.0857e18       # [cm]

# defined quantities
massComp            = massProton    # mass of the gas components
domainSize          = 10 * parsec   # cubic periodic box
ambientDens         = 10 * massComp # 10 protons per ccm 
ambientTemp         = 15            # Kelvin
adiabaticIndex      = 5.0/3.0       # mono-ionic gas (protons)
meanMolWeight       = 1             # pure Hydrogen ion

# external magnetic field flux density
Bx = 0.0
By = 0.0
Bz = 0.0

# derived quantities
polytropeConst      = boltzmannConstant * ambientTemp / massComp
ambientPres         = polytropeConst * $ambientDens
soundSpeed          = math.sqrt(polytropeConst)

maxVelocity         = machNumber * soundSpeed
turnTime            = domainSize / maxVelocity

baseForcing         = maxVelocity**2 / domainSize

maxStirTime          =  2 * turnTime
maxSimTime           = 10 * turnTime

dtinit              = turnTime / 1e6
dtmin               = turnTime / 1e10
dtmax               = turnTime / 1e2

chkPtInterval       = maxSimTime / 400


