#!/usr/bin/perl
#$number_breaks=2;
#$minimum_lines_to_skip = 200000;
chomp($ARGV[0]);
chomp($ARGV[1]);
chomp($ARGV[2]);
chomp($ARGV[3]);
$number_breaks=$ARGV[3];
$bases_per_line = -1;
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
$line=<BUFF>;
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
  my @random_numbers = map {random_int_between($min_base,$max_base)}(1..$number_breaks);#Choose a series of positions to mutate
  @random_breaks=sort {$a <=> $b} @random_numbers;
  
  my $max_lines2 = `wc -l $ARGV[1]|cut -d" " -f1`; 
  chomp $max_lines2;	
  my $max_base2 = ($max_lines2 * $bases_per_line)-($minimum_lines_to_skip*$bases_per_line);			#find out how Novel the chromosome is
  my $min_base2 = ($minimum_lines_to_skip*$bases_per_line);
  my @random_numbers2 = map {random_int_between($min_base2,$max_base2)}(1..$number_breaks);#Choose a series of positions to mutate
  @random_breaks2=sort {$a <=> $b} @random_numbers2;

###################################
###################################
#### Now add DNA breaks ###########
###################################
###################################

open (FILE , $ARGV[0])|| die "Can't read first file!\n";
$chrom1=readline(FILE); 	#skip first line while getting chromosome name
close FILE;
open (FILE , $ARGV[1])|| die "Can't read second file!\n";
$chrom2=readline(FILE); 	#skip first line while getting chromosome name
close FILE;
$chrom1=~chomp($chrom1);
$chrom1=~s/>//;
$chrom1=~ s/\t/\s/g;
$chrom1=~ s/^\s+//g;
$chrom1=~ s/\s+.+//g;
$chr1id=$chrom1;
$chr1id=~s/chr//;


$chrom2=~chomp($chrom2);
$chrom2=~s/>//;
$chrom2=~ s/\t/\s/g;
$chrom2=~ s/^\s+//g;
$chrom2=~ s/\s+.+//g;

$chr2id=$chrom2;
$chr2id=~s/chr//;

open (OUT,">>$ARGV[2]")||die "Can't write to output file!\n";
	$complete=0;

foreach $random_number(@random_breaks){
		chomp $random_number;
		my $value = $random_number/$bases_per_line;			#How many lines to skip 
		my $target_line=int(($value-0.1))+1;						#Select a line to target
		my $target_position=$bases_per_line-(($target_line*$bases_per_line)-$random_number);	 #Select poition in that line for target base
		my $temp = $target_line+1;
		$line =`awk -v stop=$temp '(NR==stop)' $ARGV[0]`;
		
		$other_chrom=shift(@random_breaks2);
		my $value2 = $other_chrom/$bases_per_line;			#How many lines to skip 
		my $target_line2=int(($value2-0.1))+1;						#Select a line to target
		my $target_position2=$bases_per_line-(($target_line2*$bases_per_line)-$other_chrom);	 #Select poition in that line for target base
		my $temp1 = $target_line2+1;
		$line2 =`awk -v stop=$temp1 '(NR==stop)' $ARGV[1]`;
		if (($line !~ /N/)&&($line2 !~ /N/)){
			chomp $line;
			my $position=$random_number;					# Adjust position
			my $ref_base=substr($line,$target_position-1,1);  # Get reference Base
			my $nextref_base=substr($line,$target_position,1);  # Get reference Base
			
			chomp $line2;
			my $position2=$other_chrom;					# Adjust position
			my $ref_base2=substr($line2,$target_position2-1,1);  # Get reference Base
			my $nextref_base2=substr($line2,$target_position2,1);  # Get reference Base
			$type=random_int_between(1,4);
			
			if($chr1id ne "" && $position  =~ m/\d+/ && $ref_base  =~ m/A|T|G|C/i && $nextref_base ne "" && $chr2id ne "" && $position2  =~ m/\d+/ && $ref_base2  =~ m/A|T|G|C/i && $nextref_base2 ne "")
			{
				if($type==1){@output=simple_translocation ($chr1id,$position,$ref_base,$nextref_base,$chr2id,$position2,$ref_base2,$nextref_base2);}
				if($type == 2){@output=single_rc_translocation($chr1id,$position,$ref_base,$nextref_base,$chr2id,$position2,$ref_base2,$nextref_base2);}
				if($type ==3){@output=double_rc_translocation($chr1id,$position,$ref_base,$nextref_base,$chr2id,$position2,$ref_base2,$nextref_base2);}
				if($type ==4){@output=reciprocal_translocation($chr1id,$position,$ref_base,$nextref_base,$chr2id,$position2,$ref_base2,$nextref_base2);}
			}
			print OUT "@output";
			$complete++;
			next}	
		}
print "Number breaks complete = $complete \n";
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
		
sub simple_translocation {
	($chr1id,$position,$ref_base,$nextref_base,$chr2id,$position2,$ref_base2,$nextref_base2)=@_;
	print "$chr1id,$position,$ref_base,$nextref_base,$chr2id,$position2,$ref_base2,$nextref_base2\n";
	$alt1=join ('',$ref_base,"[",$chr2id,":",$position2,"[");
	$alt2=join ('',"]",$chr1id,":",$position,"]",$ref_base2);
	$mate_id = join '',map $bases[rand @bases],0..5;
	$mate_id_1 =join ('',"simple_translocation_",$mate_id,"_1");
	$mate_id_2 =join ('',"simple_translocation_",$mate_id,"_2");
	return "$chr1id\t$position\t$mate_id_1\t$ref_base\t$alt1\t.\t.\tSVTYPE=BND;MATEID=$mate_id_2\n$chr2id\t$position2\t$mate_id_2\t$ref_base2\t$alt2\t.\t.\tSVTYPE=BND;MATEID=$mate_id_1\n";
	}
sub single_rc_translocation {
	($chr1id,$position,$ref_base,$nextref_base,$chr2idid,$position2,$ref_base2,$nextref_base2)=@_;
	print "$chr1id,$position,$ref_base,$nextref_base,$chr2id,$position2,$ref_base2,$nextref_base2\n";
	$alt1=join ('',$ref_base,"]",$chr2idid,":",$position2,"]");
	$alt2=join ('',$ref_base2,"]",$chr1id,":",$position,"]");
	$mate_id = join '',map $bases[rand @bases],0..5;
	$mate_id_1 =join ('',"single_rc_translocation_",$mate_id,"_1");
	$mate_id_2 =join ('',"single_rc_translocation_",$mate_id,"_2");
	return "$chr1id\t$position\t$mate_id_1\t$ref_base\t$alt1\t.\t.\tSVTYPE=BND;MATEID=$mate_id_2\n$chr2idid\t$position2\t$mate_id_2\t$ref_base2\t$alt2\t.\t.\tSVTYPE=BND;MATEID=$mate_id_1\n";
	}
sub double_rc_translocation {
	($chr1id,$position,$ref_base,$nextref_base,$chr2idid,$position2,$ref_base2,$nextref_base2)=@_;
	print "$chr1id,$position,$ref_base,$nextref_base,$chr2id,$position2,$ref_base2,$nextref_base2\n";
	$alt1=join ('',"[",$chr2idid,":",$position2,"[",$ref_base);
	$alt2=join ('',"[",$chr1id,":",$position,"[",$ref_base2);
	$mate_id = join '',map $bases[rand @bases],0..5;
	$mate_id_1 =join ('',"double_rc_translocation_",$mate_id,"_1");
	$mate_id_2 =join ('',"double_rc_translocation_",$mate_id,"_2");
	return "$chr1id\t$position\t$mate_id_1\t$ref_base\t$alt1\t.\t.\tSVTYPE=BND;MATEID=$mate_id_2\n$chr2idid\t$position2\t$mate_id_2\t$ref_base2\t$alt2\t.\t.\tSVTYPE=BND;MATEID=$mate_id_1\n";
	}
sub reciprocal_translocation {
	($chr1id,$position,$ref_base,$nextref_base,$chr2idid,$position2,$ref_base2,$nextref_base2)=@_;
	print "$chr1id,$position,$ref_base,$nextref_base,$chr2id,$position2,$ref_base2,$nextref_base2\n";
	$next_position=$position+1;
	$next_position2=$position2+1;
	$alt1=join ('',$ref_base,"[",$chr2idid,":",$position2,"[");	
	$alt2=join ('',"]",$chr1id,":",$position,"]",$ref_base2);
	$alt3=join ('',"]",$chr2idid,":",$position2,"]",$nextref_base);
	$alt4=join ('',$nextref_base2,"[",$chr1id,":",$position,"[",);
	$mate_id = join '',map $bases[rand @bases],0..5;
	$mate_id_1 =join ('',"reciprocal_translocation_",$mate_id,"_1");
	$mate_id_2 =join ('',"reciprocal_translocation_",$mate_id,"_2");
	$mate_id_3 =join ('',"reciprocal_translocation_",$mate_id,"_3");
	$mate_id_4 =join ('',"reciprocal_translocation_",$mate_id,"_4");
	return "$chr1id\t$position\t$mate_id_1\t$ref_base\t$alt1\t.\t.\tSVTYPE=BND;MATEID=$mate_id_2\n$chr2idid\t$position2\t$mate_id_2\t$ref_base2\t$alt2\t.\t.\tSVTYPE=BND;MATEID=$mate_id_1\n$chr1id\t$next_position\t$mate_id_3\t$nextref_base\t$alt3\t.\t.\tSVTYPE=BND;MATEID=$mate_id_4\n$chr2idid\t$next_position2\t$mate_id_4\t$nextref_base2\t$alt4\t.\t.\tSVTYPE=BND;MATEID=$mate_id_3\n";
	}
