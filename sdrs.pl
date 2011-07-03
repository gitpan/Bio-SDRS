#!perl -w

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

my $usage = <<EOF;
sdrs.pl [-[no]debug          ]                         input-file
        [-multiple=<float>   ] Default: $multiple
        [-ldose=<float>      ] Default: $ldose
        [-hdose=<float>      ] Default: $hdose
        [-step=<float>       ] Default: $step
        [-maxproc=<integer>  ] Default: $maxproc
        [-trim=<float>       ] Default: $trim
        [-outdir=<directory> ] Default: $outdir
EOF

if (scalar(@ARGV) != 1) {
    print STDERR $usage;
    exit(1);
}

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
if (not defined $doses) {
    die "$infile is empty.\n";
}
chomp($doses);
$doses =~ s/\s*$//;
my @doses = split (/\t/, $doses);
shift @doses;
if (scalar(@doses) == 0) {
    die "No doses specified.\n";
}

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
if (not -d $outdir) {
    my $status = system("mkdir -p $outdir");
    if ($status != 0) {
	die "Unable to create $outdir. Status = $status\n";
    }
}
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

=head1 NAME

sdrs.pl - Command line script to use Bio::SDRS to run a Signmoidal Dose Response Search.

=head1 SYNOPSIS

    sdrs.pl -multiple=1.05 \
            -step=20 \
            -ldose=0.4 \
            -hdose=25000 \
            -trim=0 \
            -outdir=results \
            data/OVCAR4_HCS_avg.txt

=head1 OPTIONS

=over
    
=item -multiple=<float>

Specifies the multiplicity factor for increasing the dose during the
search. It must be greater than one.

=item -ldose=<float>

Specifies the minimum dose for the search. Must be greater than zero.

=item -hdose=<float>

Specifies the maximum dose for the search. Must be greater than the C<ldose> above.

=item -step=<float>

This value specifies the maximum change in doses in the
search. In the search process, this module starts at the ldose
value. It tries multiplying the current dose by the C<multiple> value,
but it will only increase the dose by no more than the C<step> value
specified here.  It must be positiv.e

=item -maxproc=<integer>

Specifies the maximum number of processes to be used in the search. 

=item -trim=<float>

Trims the number of measurements used in the analysis by this factor. Must be a real number between zero and one, although zero makes no sense because no experiments would be used in the analysis.

=item -outdir=<directory>

Specifies the directory where the results are written. If this
directory does not exist, it is created.

=item -[no]debug

Controls display of debugging information. Normally should be left off.

=back

=head1 DESCRIPTION

This program provides a simple command line interface to the C<Bio::SDRS> Perl module.

=head1 INPUT FILE

There is only one input file to this script which is meant to provide
the data for one experiment. The file provides the doses of compound
used in the experiment, and it provides the measurements for each dose
for every compound used in the experiment. The file uses Tab Separated Value (TSV) format,
with each value separated from each other using the Tab character
(ASCII code 9).

The first line of the input file must contain the doses as floating
point numbers, with a comment word as the first word. For example, the
following line would describe a eight dose experiment.

Dose	0.42	1.27	3.81	11.43	34.29 	102.88	308.64	925.92	

The succeeding lines of the input file are the measured responses for
each compound tested, with no limit to the number of compounds. The
first word in each line is the name of the compound, and the remaining
lines are the responses for each dose as specified in the dose line at
the top of the file, and in the same order. Here is an example
response line consistent with the dose example above:

Aspirin	2.1	1.695	1.675	1.735	1.56	1.77	2.34	2.595

=head1 OUTPUT FILES

The output files are written to the directory specified by the C<outdir> option above. The names are composed in part from the options:

=head2 sdrs.C<multiple>.C<step>.out

The values of the sigmoidal function along with F-score for every
tested value for every compound. Each line is a tab separated value
containing the compound, dose, F-score for the fit, A, B, and D value.

=head2 sdrs.C<multiple>.C<step>.EC50.out

Provides all the EC50 data for each compound in the input file as a
tab separated record. The columns are as follows:

 1  Compound
 2  MAX         Maximum F score
 3  MIN         Minimum F score
 4  LOW         Low value for the dose.
 5  HIGH        High value for the dose.
 6  EC50        EC50 for the dose search.
 7  PVALUE      F score P value
 8  EC50RANGE   Range of dose values where F score cutoff is satisfied.
 9  PEAK        Number of peaks in the F score.
10  A           A parameter in the sigmoid function.
11  B           B parameter in the sigmoid function.
12  D           D parameter in the sigmoid function.
13  FOLD        Positive Ratio of B/A or 99999.0 if a == 0.

=head2 sdrs.sorted_probes.out

For every calculated dose in the search, list the compounds in order
of F score. Each line in this file is a TSV record where the first
field is the calculated dose, and the remaining fields are the
compound names.

=head2 sdrs.pval_FDR.out

For every calculated dose in the search, list the P Values for fit for
the compounds in order of F score exactly in the same order as the
sdrs.sorted_probes.out file.  Each line in this file is a TSV record
where the first field is the calculated dose, and the remaining fields
are the compound names. This file is meant to be used in conjunction
with the sdrs.sorted_probes.out files.

=head1 SEE ALSO

Bio::SDRS.

=head1 AUTHORS

 Ruiru Ji <ruiruji@gmail.com>
 Nathan O. Siemers <nathan.siemers@bms.com>
 Lei Ming <lei.ming@bms.com>
 Liang Schweizer <liang.schweizer@bms.com>
 Robert Bruccoleri <bruc@acm.org>

=cut
    
