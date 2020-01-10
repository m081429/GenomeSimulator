#!/usr/bin/perl
#$number_deletions=1000;
#$minimum_lines_to_skip = 200000;
chomp($ARGV[0]);
chomp($ARGV[1]);
chomp($ARGV[2]);
$number_deletions=$ARGV[2];
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
my $max_lines = ($. - 1);
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
=cut
close(BUFF);
#$minimum_lines_to_skip = 0;
### Set initial Parameters
  my @bases =("A","C","G","T");
  #my @bases =("_","_","_","_"); #For testing only
  #my $bases_per_line = 50;								#length of a fasta file
  #my $max_lines = `wc -l $ARGV[0]|cut -d" " -f1`; 
  #chomp $max_lines;	
  my $max_base = ($max_lines * $bases_per_line)-($minimum_lines_to_skip*$bases_per_line);			#find out how Novel the chromosome is
  my $min_base = ($minimum_lines_to_skip*$bases_per_line);
  my @random_numbers = map {random_int_between($min_base,$max_base)}(1..$number_deletions);#Choose a series of positions to mutate
  @random_deletions=sort {$a <=> $b} @random_numbers;
###################################
###################################
#### Now add Large Deletions ######
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
	my $min_deletion_size=1000;
	my $max_deletion_size=10000;

foreach $random_number(@random_deletions){
		chomp $random_number;
		my $value = $random_number/$bases_per_line;			#How many lines to skip 
		my $target_line=int(($value-0.1))+1;						#Select a line to target
		my $target_position=$bases_per_line-(($target_line*$bases_per_line)-$random_number);	 #Select poition in that line for target base
		my $temp = $target_line+1;
		$line =`awk -v stop=$temp '(NR==stop)' $ARGV[0]`;
			if ($line !~ /N/){
			chomp $line;
			$deletion_length = random_int_between($min_deletion_size,$max_deletion_size);
			my $position=$random_number;					# Adjust position
			my $ref_base=substr($line,$target_position-1,1);  # Get reference Base
			my $deletion_length=$deletion_length+1;
			my $end=$deletion_length+$position;
			if($ref_base eq ""){next};
			if($position =~ m/\d+/ && $ref_base =~ m/A|T|G|C/i)
			{
			print OUT "$chrid\t$position\t.\t$ref_base\t<DEL>\t.\t.\tSVTYPE=DEL;END=$end;SVLEN=$deletion_length\n";
			$complete++;
			}
			next}	
		}
print"Number Novel deletions complete = $complete \n";
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
