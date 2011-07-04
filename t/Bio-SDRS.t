# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Bio-SDRS.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More qw(no_plan);
BEGIN { use_ok('Bio::SDRS') };

my $multiple = 1.05;
my $ldose = 0.4;
my $hdose = 25000;
my $step = 20;
my $maxproc = 4;
my $trim = 0;
my $outdir = ".";
my $debug = 0;

my $sdrs = Bio::SDRS->new();
ok($sdrs, "Object created.");
ok($sdrs->multiple($multiple) == $multiple, "Multiple set");
ok($sdrs->ldose($ldose) == $ldose, "ldose set");
ok($sdrs->hdose($hdose) == $hdose, "hdose set");
ok($sdrs->step($step) == $step, "step set");
ok($sdrs->maxproc($maxproc) == $maxproc, "maxproc set");
ok($sdrs->trim($trim) == $trim, "trim set");
ok($sdrs->debug($debug) == $debug, "debug set");

my $infile = "t/OVCAR4_HCS_avg.txt";

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
my $file = "t/sdrs.$multiple.$step.out";
open (OUT, ">$file") ||
    die "Unable to open $file: $!\n";
print OUT $sdrs->scandata;
close OUT;
open (OUT, ">t/sdrs.$multiple.$step.EC50.out") ||
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

open (SORT, ">t/sdrs.sorted_probes.out") ||
    die "can not open sorted probeset output file: $!\n";
open (PVAL, ">t/sdrs.pval_FDR.out") ||
    die "can not open p value output file: $!\n";

foreach my $dose ($sdrs->score_doses) {
    my $dose_st = sprintf("%.5f", $dose);
    print SORT "${dose_st}\t", join("\t", $sdrs->sorted_assays_by_dose($dose)), "\n";
    print PVAL "${dose_st}\t", join("\t", $sdrs->pvalues_by_dose($dose)), "\n";
}

close SORT;
close PVAL;
$ENV{"LC_ALL"} = "C";
foreach my $f (('sdrs.1.05.20.EC50.out',
		'sdrs.1.05.20.out',
		'sdrs.pval_FDR.out',
		'sdrs.sorted_probes.out')) {
    ok(system("sort t/${f} >t/${f}.srt") == 0, "Sort ${f}");
    &compare_files("t/${f}.srt", "t/ref.${f}.srt");
}

sub compare_files {
    my $file1 = shift;
    my $file2 = shift;
    my $limit = 1000;

    local (*IN1, *IN2);

    ok(open (IN1, "<$file1"), "Open $file1");
    ok(open (IN2, "<$file2"), "Open $file2");
    my $count = 0;
    while (my $line1 = <IN1>) {
	my $line2 = <IN2>;
	if (not defined $line2) {
	    ok(0, "$file1 bigger than $file2");
	    last;
	}
	$count++;
	if ($line1 eq $line2) {
	    ok(1, "Lines $count match");
	}
	else {
	    if ($limit-- > 0) {
		ok($line1 eq $line2,
		   "line $count match: line1 = $line1  line2 = $line2");
	    }
	}
    }
    my $line2 = <IN2>;
    if (defined $line2) {
	ok(0, "$file2 bigger than $file1");
    }
    close IN1;
    close IN2;
}


