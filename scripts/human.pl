#!/usr/bin/perl -w

#
# Build an Ebwt out of a human genome build in the specified directory
#

use Getopt::Std;

# Trim whitespace from a string argument
sub trim($) {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

# Print usage message
sub usage() {
	print "Usage: perl human.pl [-b <bmax>] [-c <chr_list>] [-h <genome_dir>] <run_name> [<outdir>]\n";
	print "Runs ebwt_build on some or all of the human genome\n";
	print "   -b <int_list>       set bucket size to each elt of list for blockwise SA builder\n";
	print "   -c <chr_list>       comma-delimited list of chromosomes to include (e.g. 3,4,7)\n";
	print "   -h <genome_dir>     set base dir for human genome (overrides HS_GENOME_DIR env var)\n";
	print "   -d                  dry-run; read command-line arguments but do not execute command\n";
	print "   -a                  run with-asserts versions of programs\n";
	print "   -p                  run packed-string versions programs\n";
	print "   -r \"<args>\"         pass <args> to ebwt_build\n";
	exit 1;
}

my %options=();
getopts("c:b:h:dvpar:",\%options);

my $verbose = 0; $verbose = $options{v} if defined $options{v};
my $dry_run = 0; $dry_run = $options{d} if defined $options{d};
my $packed  = 0; $packed  = $options{p} if defined $options{p};

# Set bmax, possibly from command line
my @bmax = (1024*1024);
@bmax = split(/,/, $options{b}) if defined $options{b};

# Set human genome dir; possibly from environment variable, possibly
# from command line
my $hs_dir = "";
$hs_dir = $ENV{HS_GENOME_DIR} if defined $ENV{HS_GENOME_DIR};
$hs_dir = $options{h} if defined $options{h};

my $asserts = "";
$asserts = "-with-asserts" if defined $options{a};

my $args = "";
$args = $options{r} if defined $options{r};

# Set chromosomes, possibly from command line
my @chrs = ( "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",
             "9", "10", "11", "12", "13", "14", "15", "16",
            "17", "18", "19", "20", "21", "22",  "X",  "Y");
@chrs = split(/,/, $options{c}) if defined $options{c};

# Check that all required arguments are specified
if($hs_dir eq "") {
	print "No human-genome base directory specified\n"; usage();
}
if($#chrs == -1) {
	print "No chromosomes specified\n"; usage();
}
if($#ARGV < 0) {
	print "Missing arguments\n"; usage();
}

# Print 
print "Using bmax: @bmax\n" if $verbose;
print "Using chromosomes: @chrs\n" if $verbose;
print "Using human genome dir: $hs_dir\n" if $verbose;
print "Using ebwt_build arguments: $args\n" if $verbose;

my $runname = $ARGV[0];
if($runname eq "") {
	print "Missing runname\n"; usage();
}
print "Using run name: $runname\n" if $verbose;

my $outdir = ".";
$outdir = $ARGV[1] if defined($ARGV[1]);
if($runname eq "") {
	print "Missing outdir\n"; usage();
}
print "Using output directory: $outdir\n" if $verbose;

print "Asserts enabled\n" if $asserts;

# Build a text file list that includes every chromosome
# TODO: Mitochondrial?
my $inputs = '';
foreach my $idx (@chrs) {
    $inputs .= "$hs_dir/hs_ref_chr$idx.mfa,";
}
chop($inputs); # Remove trailing comma

print "Inputs: $inputs\n" if $verbose;

my $profile = 1;
my $tee = 1;

if($packed) {
	$packed = "_packed";
} else {
	$packed = "";
}

# Execute command for each specified bmax
foreach my $b (@bmax) {
	# Execute ebwt_build
	my $cmd = "";
	$cmd .= "perl scripts/profile_wrap.pl $runname.$b $outdir " if $profile;
	$cmd .= "./ebwt_build$packed$asserts --profile -v --bmax $b $args $inputs $outdir/$runname";
	# Don't include bucket name in output file for now; we'd prefer to
	# overwrite previous .ebwts just for space reasons
	#$cmd .= "./ebwt_build$packed$asserts --profile -v --bmax $b $args $inputs $outdir/$runname.$b";
	$cmd .= " | tee $outdir/build.$runname.$b.out 2>&1" if $tee;
	print "$cmd\n";
	if(!$dry_run) {
		my $exitlevel = system("$cmd");
		if($exitlevel != 0) {
			print "Non-zero exitlevel from ebwt_build: $exitlevel\n";
			exit $exitlevel;
		}
	}
}
