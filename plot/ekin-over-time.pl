#!/usr/bin/env perl

use v5.10;
use strict;
use warnings;

my $machnumber  = shift @ARGV;
my $turntime    = shift @ARGV;
my @datafiles   = ();
my @legends     = ();

while (@ARGV) {
    push @datafiles, shift @ARGV;
    push @legends, shift @ARGV;
}

my $cmds = <<"END";
#set terminal svg enhanced dashed
set terminal pngcairo enhanced dashed crop
set encoding utf8

set size ratio 0.6

mach		= $machnumber
crosstime	= $turntime 
datafiles	= \"@datafiles\"
legends		= \"@legends\"

scale_x(x) = x / crosstime
scale_y(x) = x

set grid

set xrange[0:5]

set xtics 0.5

set title "Total Kinetic Energy {/:Italic 𝓚} over crossing time {/:Italic t / t_c}"

set xlabel 'crossing time {/:Italic t / t_c}   →'
set ylabel 'kinetic energy {/:Italic 𝓚}   →'

plot for [i=1:words(datafiles)] word(datafiles,i) \\
		u (scale_x(\$3)):(scale_y(\$4)) w l lw 2 t word(legends,i)
END

open GP, "| gnuplot";
print GP $cmds;
