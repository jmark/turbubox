#!/usr/bin/env perl

use v5.10;
use strict;
use warnings;

use Config::Tiny;
use Getopt::Long qw/GetOptionsFromArray/;

my %config = (
    energy_policy_tag         => 'devel_markert',
    minimize_time_to_solution => 'yes',
    
    notify_user      => 'markert@ph1.uni-koeln.de',
    notification     => 'always',
    
    class            => 'micro',
    job_type         => 'parallel',
    output           => 'plot.log',
    error            => 'plot.log',
    
    job_name         => 'micro test',
    initialdir       => $ENV{PWD},
    
    island_count     => 1,
    node             => 1,
    total_tasks      => 1,
    wall_clock_limit => '01:30:00',
    shell            => '/bin/bash',
    node_usage       => 'not_shared',
);

my @argv = ();
while (my $el = shift @ARGV) {
    last if $el eq '::';
    push @argv, $el;
}

my @opts = ();
while (my $el = shift @ARGV) {
    last if $el eq '::';
    push @opts, $el;
}

my @cmd = @ARGV;

my @files = ();
my $output = '';

GetOptionsFromArray(\@argv, 'file=s' => \@files, 'out=s'  => \$output);

for my $file (@files) {
    my $ini = Config::Tiny->read($file);
    %config = (%config, %{$ini->{_}});
}

GetOptionsFromArray(\@opts, map {;"$_=s" => \$config{$_};} keys %config);

say '#!/bin/bash';
say '';

while (my ($k,$v) = each %config) {
    $k =~ s/_DOT_/./g; # Getopt::Long has problems with '.' in keys
    say "#@ $k = $v";
}

say '#@ queue';
say '';
say 'source /etc/profile';
say 'source /etc/profile.d/modules.sh';
say '';
say join ' ', ('eval', map {;"'$_'";} @cmd);
