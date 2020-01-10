#!/usr/bin/perl
#### Genome_mutator.pl ###
## hart.steven@mayo.edu
##modified:naresh.prodduturi@mayo.edu
use Getopt::Long;
use Benchmark;

my $inputChr = '';
my $outputVCF = 'snp_mutator_output.vcf';
my $mutations = 50000;
my $insertions = '';
my $deletions = '';
GetOptions(
	'f=s' => \$inputChr,
	'o:s' => \$outputVCF,
	'm:i' => \$mutations,
	'i:i' => \$insertions,
	'd:i' => \$deletions)|| die usage();
sub usage {
	print "\n\tusage: SNP_mutator.pl -f <Chr##.fa> \n";
	print "\t-o\tOutput File. Default = snp_mutator_output.vcf\n";
	print "\t-m\tTotal number of Desired Mutations\n";
	print "\t-i\tTotal number of Desired Insertions\n";
	print "\t-d\tTotal number of Desired Deletions\n\n";
	}
	
if($inputChr eq ''){die "No Chr File Provided! \n"}
print "Output: $outputVCF\n\n";

my @bases =("A","C","G","T");
#die  "$bases_per_line\n";
if($insertions == ''){ $insertions=int(($mutations-0.1)/10)+1; }
if($deletions == ''){ $deletions=int(($insertions-0.1)/2)+1; }

open (OUT, ">$outputVCF") || die "Can't write to output to $outputVCF!\n";
#print OUT "all\n";
	
###################################
#### Start by Introducing SNVs ####
###################################
$random_places=($mutations+$insertions+$deletions);

## If for whatever reason, numbers are used in place of ChrX,Y,M --> replaced
#$inputChr=~s/chr23(.*)chr24(.*)chr25/chrX${1}chrY${1}chrM/;
## Find where the Genetic Code begins, avoid N's
$minimum_lines_to_skip =0;
$t0 = Benchmark->new;  ### ---------> Benchmark
open(INPUT,"<$inputChr") || die "Unable to Open $inputChr\n";
$bases_per_line =-1;
while(<INPUT>)
{
	next if($_ =~ /^\>/); ##Skip headers
	if($bases_per_line != -1 && $bases_per_line != length($_)-1)
	{
		die "unequal fasta lines length\n";
	}
	else
	{
		$bases_per_line=length($_)-1;
	}
	if($_ !~ m/A|C|G|T/i){$minimum_lines_to_skip++;}
	else{last;}
}
while(<INPUT>)
{}
my $max_lines = ($. - 1); ##Conitue to Last Line subtract header

#print "Working Lines: $minimum_lines_to_skip - $max_lines\n";
$t1 = Benchmark->new;
$td = timediff($t1, $t0);
print "Read Chr File took: ",timestr($td),"\n";

### Set initial Parameters #######################################
my $max_base = ($max_lines * $bases_per_line);	#find out how large the chromosome is
my $min_base = ($minimum_lines_to_skip*$bases_per_line);
print "Max=$max_base\tMin=$min_base\t$minimum_lines_to_skip\t$max_lines\n";
my @random_numbers=0;
my @random_numbers = map {random_int_between($min_base,$max_base)}(1..$random_places);#Choose a series of positions to mutate
#print "check 1 @random_numbers\n";
my @random_snps=@random_numbers[0..$mutations-1];
my @random_insertions=@random_numbers[$mutations..$mutations+$insertions-1];
my @random_deletions=@random_numbers[$mutations+$insertions..$random_places-1];
#die "check2 @random_snps\n @random_insertions\n @random_deletions\n"; 
  @random_snps=sort {$a <=> $b} @random_snps;
  @random_insertions=sort {$a <=> $b} @random_insertions;
  @random_deletions=sort {$a <=> $b} @random_deletions;
#  print ">>>>>>>@random_snps\n";
#################################################################

$t2 = Benchmark->new;
$td2 = timediff($t2, $t1);
print "Generate Random Arrays took: ",timestr($td2),"\n\n";

seek INPUT,0,0;  #Reset to file beginning without Open/Close.
$chrom = <INPUT>;
my $chrid;
if($chrom =~ /chr(\w+)/i){$chrid=$1;}
#die "check3 $chrid\n";
###################################
#### Create SNV positions #########
###################################
$complete=0; my $n=0; $.=0; ## reset line counter.
$target_line = lineOfSNP( $random_snps[$n] );
#$last_line = lineOfSNP( $random_snps[-1] )+1;
#die "$target_line\t$last_line\n";
#print "positions @random_snps\n";
while(<INPUT>){
	if($. == $target_line ){
		my $inLinePos= ($bases_per_line - (($target_line*$bases_per_line)-$random_snps[$n]))-1;
		#die "$random_snps[$n] $inLinePos\n";
		@nucleo = split('', $_);
		if( !defined($nucleo[$inLinePos]) )
		{ ##just in case
			#print "SNP REJECTED : Not Defined POS:$random_snps[$n] LINE:$target_line Pos:$inLinePos\n";
			$n++; 
			$target_line = lineOfSNP( $random_snps[$n] );
			redo;
		} 
		if( $nucleo[$inLinePos] =~ /N/ || $nucleo[$inLinePos] =~ /\n/ )
		{ ### Skip if base is N or nextLine
			#print "SNP REJECTED : N or \\n POS:$random_snps[$n] LINE:$target_line\n";
			$n++; 
			$target_line = lineOfSNP( $random_snps[$n] );
			redo;
		} 
		$base=randNonSelfBase($nucleo[$inLinePos]);	
		if($nucleo[$inLinePos] eq $base)
		{
			redo;
		}
		print OUT "$chrid\t$random_snps[$n]\t.\t$nucleo[$inLinePos]\t$base\t.\t.\t.\n";
		#print OUT "$chrid\t".$.."=$target_line:$target_position\t$random_snps[$n]\t.\t$nucleo[$inLinePos]\t$base\t.\t.\t.\n";
		$complete++; 
		$n++; 
		$target_line = lineOfSNP( $random_snps[$n] );
		if($. == $target_line){redo;}
	}
}
print"SNVs requested: ".@random_snps."\tcompleted: $complete\n";

$t3 = Benchmark->new;
$td3 = timediff($t3, $t2);
print "SNV positions took: ",timestr($td3),"\n\n\n\n\n\n";
#die;
###################################
#### Now add insertions ###########
###################################
seek INPUT,0,0;  #Reset to file beginning without Open/Close.
$header = <INPUT>;
$complete=0; my $n=0;
$.=0; ## reset line counter.
$target_line = lineOfSNP( $random_insertions[$n] );
$last_line = lineOfSNP( $random_insertions[-1] )+1;
#print "insertion positions @random_insertions\n";
while(<INPUT>){
	if($. == $target_line){
		my $inLinePos= ($bases_per_line - (($target_line*$bases_per_line)-$random_insertions[$n]))-1;
		@nucleo = split('', $_);
		if( !defined($nucleo[$inLinePos]) ){ ##just in case
			#print "INSERTION REJECTED: Not Defined POS:$random_insertions[$n] LINE:$target_line SNP REJECTED Pos:$inLinePos\n";
			$n++; $target_line = lineOfSNP( $random_insertions[$n] );
			redo;
		} 
		elsif( $nucleo[$inLinePos] =~ /N/ || $nucleo[$inLinePos] =~ /\n/){ ### Skip if base is N or nextLine
		#print "INSERTION REJECTED : N or \n POS:$random_insertions[$n] LINE:$target_line\n";
			$n++; $target_line = lineOfSNP( $random_insertions[$n] );
			redo;
		} 
		$insertion=randInsertionArray($nucleo[$inLinePos]);
		if($insertion !~ m/^$nucleo[$inLinePos]/)
		{
			redo;
		}	
		print OUT "$chrid\t$random_insertions[$n]\t.\t$nucleo[$inLinePos]\t$insertion\t.\t.\t.\n";
		$complete++; 
		$n++; 
		$target_line = lineOfSNP( $random_insertions[$n] );
		if($. == $target_line){redo;}
	}
}
print"Inserts requested: ".@random_insertions."\tcompleted: $complete\n";

$t4 = Benchmark->new;
$td4 = timediff($t4, $t3);
print "Insertions took: ",timestr($td4),"\n\n\n\n\n\n";
#die;
###################################
#### Now add deletions ############
###################################
$complete=0;my $n=0;
seek INPUT,0,0;  #Reset to file beginning without Open/Close.
$header = <INPUT>;
$.=0; ## reset line counter.
$target_line = lineOfSNP( $random_deletions[$n] );
$last_line = lineOfSNP( $random_deletions[-1] )+1;

while(<INPUT>){
	if($. == $target_line){
		my $inLinePos= ($bases_per_line - (($target_line*$bases_per_line)-$random_deletions[$n]))-1;
		@nucleo = split('', $_);
		if( !defined($nucleo[$inLinePos]) ){ ##just in case
			#print "DELETION REJECTED : Not Defined POS:$random_deletions[$n] LINE:$target_line\n";
			$n++; $target_line = lineOfSNP( $random_deletions[$n] );
			redo;
		} 
		elsif( $nucleo[$inLinePos] =~ /N/ || $nucleo[$inLinePos] =~ /\n/){ ### Skip if base is N or nextLine
			#print "DELETION REJECTED: N or \\n POS:$random_deletions[$n] LINE:$target_line\n";
			$n++; $target_line = lineOfSNP( $random_deletions[$n] );
			redo;
		} 
		$deletionLn=randDeletionArray($inLinePos);
		$deletion  = substr $_, $inLinePos, $deletionLn;
		$deletion=~ s/\n//g;
		#if($deletion !~ m/\w/)
		#{
			#print "check main $inLinePos\n";
			
		#}
		if($deletion =~ m/^\w/ && $nucleo[$inLinePos] =~ m/^\w/)
		{
		print OUT "$chrid\t$random_deletions[$n]\t.\t$deletion\t$nucleo[$inLinePos]\t.\t.\t.\n";# edited by SNH to comply with VCF spec
	#	print OUT "$chrid\t$random_deletions[$n]\t.\t$deletion\t$nucleo[$inLinePos]\t.\t.\t.\n";
		$complete++; $n++; 
		$target_line = lineOfSNP( $random_deletions[$n] );
		}
		if($. == $target_line){redo;}
	}
}
print"Deletes requested: ".@random_deletions."\tcomplete: $complete \n";

$t5 = Benchmark->new;
$td5 = timediff($t5, $t4);
print "Create Deletions took: ",timestr($td5),"\n\n";

close OUT;


########################################################################################################################
### SUBROUTINES
########################################################################################################################
sub random_int_between {
		my($min, $max) = @_;
		# Assumes that the two arguments are integers themselves!
		return $min if $min == $max;
		($min, $max) = ($max, $min)  if  $min > $max;
		return $min + int rand(1 + $max - $min);
		}
		
sub lineOfSNP {
	my $value = $_[0]/$bases_per_line;
	my $target_line=int($value)+1;
	#print "\t POS: $_[0] -> Line: $target_line\n";
	return $target_line;
}

sub randNonSelfBase{
	my @remaining;
	if(uc($_[0]) eq "A"){ @remaining =("C","G","T"); }
	elsif(uc($_[0]) eq "T"){ @remaining =("C","G","A"); }
	elsif(uc($_[0]) eq "C"){ @remaining =("A","G","T"); }
	elsif(uc($_[0]) eq "G"){ @remaining =("C","A","T"); }
	$b = $remaining[rand @remaining];
	return $b;
}

sub randInsertionArray{
	my $max_base = $bases_per_line;
	my $min_base = 1;
	@insertion_length1 = map {random_int_between($min_base,$max_base)}(1..2);
	my $max_base = 5;
	my $min_base = 1;
	@insertion_length2 = map {random_int_between($min_base,$max_base)}(1..60);	# Make small insertions 30x more likely than large ones
	@insertion_length3 =(1) x 120;
	@insertion_length=(@insertion_length1,@insertion_length2,@insertion_length3);
	my $insertion_length=$insertion_length[rand @insertion_length];
	my $insertion = join '',map $bases[rand @bases],0..$insertion_length;
	return $_[0]."".$insertion;
}

sub randDeletionArray{
	my $max_base = $bases_per_line-$_[0];
	my $min_base = 1;
	@deletion_length1 = map {random_int_between($min_base,$max_base)}(1..5);
	@deletion_length3 =(1) x 5;
	@deletion_length=(@deletion_length1,@deletion_length3);
	my $deletion_length=$deletion_length[rand @deletion_length]+1;
	return 	$deletion_length;
}


