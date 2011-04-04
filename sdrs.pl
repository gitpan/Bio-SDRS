#!/usr/bin/perl -w

use strict;
use lib "lib";
use Getopt::Long;
use Bio::SDRS;

my $multiple = 1.18;
my $ldose = 0.17;
my $hdose = 30000;
my $step = 60;
my $maxproc = 20;
my $trim = 0.6;
my $outdir = ".";
my $debug = 0;
&GetOptions("multiple=f" => \$multiple,
	    "ldose=f" => \$ldose,
	    "hdose=f" => \$hdose,
	    "step=f" => \$step,
	    "maxproc=i" => \$maxproc,
	    "trim=f" => \$trim,
	    "outdir=s" => \$outdir,
	    "debug!" => \$debug);
my $sdrs = Bio::SDRS->new();
$sdrs->multiple($multiple);
$sdrs->ldose($ldose);
$sdrs->hdose($hdose);
$sdrs->step($step);
$sdrs->maxproc($maxproc);
$sdrs->trim($trim);
$sdrs->debug($debug);

my $infile = $ARGV[0];

open (IN, $infile) || die "can not open infile $infile: $!\n";
my $doses = <IN>;
chomp($doses);
$doses =~ s/\s*$//;
my @doses = split (/\t/, $doses);
shift @doses;

$sdrs->doses(@doses);
my $count = 0;
while (<IN>) {
    chomp;
    $count++;
    my ($assay, @data) = split (/\t/, $_);
    $sdrs->set_assay($assay, @data);
}
close IN;
$sdrs->calculate;
my $file = "$outdir/sdrs.$multiple.$step.out";
open (OUT, ">$file") ||
    die "Unable to open $file: $!\n";
print OUT $sdrs->scandata;
close OUT;
open (OUT, ">$outdir/sdrs.$multiple.$step.EC50.out") ||
    die "can not open EC50 output file: $!\n";
foreach my $assay ($sdrs->assays) {
    print OUT "$assay";
    foreach my $prop (('MAX', 'MIN', 'LOW', 'HIGH', 'EC50',
		       'PVALUE', 'EC50RANGE', 'PEAK', 'A', 'B',
		       'D', 'FOLD')) {
	print OUT "\t", $sdrs->ec50data($assay, $prop);
    }
    print OUT "\n";
}
close OUT;

open (SORT, ">$outdir/sdrs.sorted_probes.out") ||
    die "can not open sorted probeset output file: $!\n";
open (PVAL, ">$outdir/sdrs.pval_FDR.out") ||
    die "can not open p value output file: $!\n";

foreach my $dose ($sdrs->score_doses) {
    print SORT "$dose\t", join("\t", $sdrs->sorted_assays_by_dose($dose)), "\n";
    print PVAL "$dose\t", join("\t", $sdrs->pvalues_by_dose($dose)), "\n";
}

close SORT;
close PVAL;



