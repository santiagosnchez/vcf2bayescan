use List::MoreUtils uniq;
my $gzip = 0;
my $vcf;
my $popfile;
my $al;
my %pops=();
my @exclude=();
my @indnms=();
my %data=();
my $usage = "Usage:
perl vcf2bayescan.pl -p <popfile> -v <vcf-file> -e <exclude1:[exclude2]>

example:
perl vcf2bayescan.pl -p pops.txt -v snps.vcf -e indv10:indiv12

pop file example:
ind1<tab>pop1
ind2<tab>pop1
ind3<tab>pop2
ind4<tab>pop2
ind5<tab>pop3
ind6<tab>pop3\n";

# arguments

if ((scalar(@ARGV) == 0) or (grep { /-+he{0,1}l{0,1}p{0,1}/ } @ARGV)){
	die $usage;
}

if (my ($indV) = grep { $ARGV[$_] =~ /^-v$/ } 0 .. $#ARGV){
	$vcf = $ARGV[$indV+1];
	if ($vcf =~ m/\.gz$/){
		print "VCF is compressed\n";
		open(VCF, "gunzip -c $vcf |") or die "can't open gunzip pipe to $vcf\n";
	} else {
		open(VCF, "<", $vcf) or die "can't open $vcf\n";
	}
} else {
	die $usage;
}

if (my ($indP) = grep { $ARGV[$_] =~ /^-p$/ } 0 .. $#ARGV){
	$popfile = $ARGV[$indP+1];
	open(POP, "<", $popfile) or die "can't open $popfile\n";
} else {
	die $usage;
}

if (my ($indE) = grep { $ARGV[$_] =~ /^-e$/ } 0 .. $#ARGV){
	if ($ARGV[$indE+1] =~ m/:/){
		@exclude = split(/:/, $ARGV[$indE+1]);
	} else {
		push @exclude, $ARGV[$indE+1];
	}
}

# read and process pop file

while (<POP>){
	chomp $_;
	@line = split /\t/, $_;
	$pops{$line[0]} = $line[1];
}
close POP;

# Get populations and sort them

my @upops = sort { $a cmp $b } uniq ( values %pops );
print "found ", scalar @upops, " populations\n";

# read and process VCF file

local $| = 1;
while(<VCF>){
	if (/^##/){
		$head .= $_;
	}
	elsif (/^#CHROM/){
		chomp($_);
		@HEAD = split /\t/, $_;
		($format_ind) = grep { @HEAD[$_] =~ m/FORMAT/ } 0 .. $#HEAD;
		@indnms = @HEAD[ ($format_ind+1) .. $#HEAD ];
		@indidx = 0 .. $#indnms;
		if (scalar(@exclude) != 0){
			@exidx = map { $ex = $_;
				       ($x) = grep { $indnms[$_] =~ /$ex/ } 0 .. $#indnms;
				       splice(@indidx,$x,1);
				       splice(@indnms,$x,1) } @exclude;
		}
		$cline = $. + 1;
		$al = 0;
	}
	else {
		++$al;
		chomp($_);
		@DAT = split /\t/, $_;
		for $pop (@upops){
			@get_ind_pop=();
			map { push @get_ind_pop, $_ if ($pops{$_} eq $pop) } keys %pops;
			@myind = getind(\@HEAD,\@get_ind_pop);
			@site = @DAT[@myind];
			map { s/:.*// } @site;
			map { s/\/|\|//g } @site;
			$gentmp = join("", @site);
			$gentmp =~ s/\.//g;
			%COUNTS = count(split //, $gentmp);
			if (!$COUNTS{0}){ $COUNTS{0} = 0 };
			if (!$COUNTS{1}){ $COUNTS{1} = 0 };
			$total = eval(join("+",values %COUNTS));
			$numal = scalar(values %COUNTS);
			$outline = sprintf("% 5s % 5s % 5s % 5s % 5s\n", $al, $total, $numal, $COUNTS{0}, $COUNTS{1});
			$data{$pop} .= $outline;
		}
		print "processing SNP $al\r";
	}
}
close VCF;
local $| = 0;

print "\nsaving to file: bayescan.txt\n";

open(OUT, "> bayescan.txt");
print OUT "[loci]=$al\n\n";
print OUT "[populations]=",scalar(@upops),"\n\n";
foreach(@upops){
	print OUT "[pop]=$_\n";
	print OUT $data{$_};
}
close OUT;
print "done. :-)\n";

# subroutines

sub getind {
	my ($line,$names) = @_;
	my @ind = map { $el = $_; grep { $$line[$_] =~ m/^$el$/ } 0 .. $#$line } @$names;
	return @ind;
}

sub count {
	my @arr = @_;
	my %co;
	++$co{$_} for @arr;
	return %co;
}
