#!/usr/bin/perl
#### Genome_mutator.pl ###
#$mutations=1000;	#number of mutations to make
#$insertions=int(($mutations-0.1)/10)+1;
#$deletions=int(($insertions-0.1)/2)+1;
$number_Novel_insertions=0;
#$random_places=($mutations+$insertions+$deletions);
#$minimum_lines_to_skip = 200000;
chomp($ARGV[0]);
chomp($ARGV[1]);
chomp($ARGV[2]);
$number_Novel_insertions= $ARGV[2];
my $bases_per_line = -1;
if($ARGV[0] =~ m/chr23/)
{
        $ARGV[0] =~ s/chr23/chrX/g;
}
if($ARGV[0] =~ m/chr24/)
{
        $ARGV[0] =~ s/chr24/chrY/g;
}
if($ARGV[0] =~ m/chr26/)
{
        $ARGV[0] =~ s/chr26/chrM/g;
}

$minimum_lines_to_skip =0;
open(BUFF,"$ARGV[0]") or die "no file found\n";
$line=<BUFF>;
=head
$k =0;
while($k == 0)
{
        $line=<BUFF>;
        if($line =~ m/^N/)
        {
                $minimum_lines_to_skip++;
        }
        else
        {
                $k = 1;
        }

}
close(BUFF);
=cut
while(<BUFF>)
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
while(<BUFF>){}
my $max_lines = ($. - 1); ##Conitue to Last Line subtract header
$tempdir = "INS_DIR";
#creating temp directory if not exists
mkdir "$tempdir", unless -d "$tempdir";

#$minimum_lines_to_skip = 0;
### Set initial Parameters
  my @bases =("A","C","G","T");
  #my @bases =("_","_","_","_"); #For testing only
  #my $bases_per_line = 50;								#length of a fasta file
  #my $max_lines = `wc -l $ARGV[0]|cut -d" " -f1`; 
  #chomp $max_lines;	
  my $max_base = ($max_lines * $bases_per_line)-($minimum_lines_to_skip*$bases_per_line);			#find out how Novel the chromosome is
  my $min_base = ($minimum_lines_to_skip*$bases_per_line);
  my @random_numbers = map {random_int_between($min_base,$max_base)}(1..$number_Novel_insertions);#Choose a series of positions to mutate
  @random_insertions=sort {$a <=> $b} @random_numbers;
###################################
###################################
#### Now add Novel Insertions #####
###################################
###################################

open (FILE , $ARGV[0])|| die "Can't read FASTA file!\n";
$chrom=readline(FILE); 	#skip first line while getting chromosome name
close FILE;
$chrom=~chomp($chrom);
$chrom=~s/>//;
$chrom=~ s/\t/\s/g;
$chrom=~ s/^\s+//g;
$chrom=~ s/\s+.+//g;
$chrom=~ s/chr//g;
$chrid=$chrom;
$chrid=~s/chr//;


open (OUT,">>$ARGV[1]")||die "Can't write to output file!\n";
	$complete=0;
	my $min_insertion_size=1000;
	my $max_insertion_size=10000;

foreach $random_number(@random_insertions){
		chomp $random_number;
		my $value = $random_number/$bases_per_line;			#How many lines to skip 
		my $target_line=int(($value-0.1))+1;						#Select a line to target
		my $target_position=$bases_per_line-(($target_line*$bases_per_line)-$random_number);	 #Select poition in that line for target base
		my $temp = $target_line+1;
		$line =`awk -v stop=$temp '(NR==stop)' $ARGV[0]`;
			if ($line !~ /N/){
			chomp $line;
			$insertion_length = random_int_between($min_insertion_size,$max_insertion_size);
			$insertion_sequence = join '',map $bases[rand @bases],0..$insertion_length;
			my $position=$random_number;					# Adjust position
			my $ref_base=substr($line,$target_position-1,1);  # Get reference Base
			my $insertion_length=$insertion_length+1;
			#$ref_base =~ s/[\t,\n,\t, ]//g;
			if($position =~ m/\d+/ && $ref_base =~ m/A|T|G|C/i)
			{
			print  OUT "$chrid\t$position\tINS.$complete\t$ref_base\t<$chrom.INS.$complete.fa>\t\.\t\.\tSVTYPE=INS;SVLEN=$insertion_length\n";
			open (INSERT,">$tempdir/$chrom.INS.$complete.fa")||die "Can't write to output file!\n";
			print INSERT ">$chrom.INS.$complete\n$insertion_sequence\n";
			$complete++;
			}
			next}	
		}
print"Number Novel Insertions complete = $complete \n";
close OUT;
close INSERT;
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
