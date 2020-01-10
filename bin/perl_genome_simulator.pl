#######################################
#Naresh Prodduturi
#10-05-2011
#Genome Simulator
#
#######################################
#reading input param & files
#!/usr/bin/perl
=head
use Cwd 'abs_path';
my $line = abs_path($0);
chomp $line;
my @DR_array = split('/',$line);
pop(@DR_array);
#pop(@DR_array);
my $dir = join("/",@DR_array);
#$dir = "$dir/Modules/";
#print "$dir/Modules/\n";
#die;
#use lib "/data1/bsi/BORA_processing/devel/genome_simulator/modifiedscripts/Modules";
#use lib "$dir/Modules";
#use vars qw($Bin $Script);
#BEGIN { ($Bin, $Script) = split /([^\/\\]+)$/, $0 }
#use lib $Bin/Modules;
=cut
use FindBin;
use lib "$FindBin::Bin/../Modules";
use File::ReadBackwards; 
use warnings;
use Cwd 'abs_path';

#use strict;
my $refgenome_dir;
my $input_vcf_like;
my $ins_files_dir;
my $vcf_chr;
my $vcf_pos;
my $vcf_id;
my $vcf_ref;
my $vcf_alt;
my $vcf_qual;
my $vcf_filter;
my $vcf_info;
my $svtype;
my $ins_filename;
my $linenum;
my $svlen;
my $temp;
my $file;
my $i;
my $chr_end_cord;
my $buff;
my $line;
my $length;
my $maternal;
my $paternal;
my $tempdir;
my $index;
my $flag;
my $temp2;
my $temp3;
my $eventstart;
my $eventstop;
my $cnv;
my $maincnv;
my $maincnvnum;
my $ref1;
my $ref2;
my $flag2;
my $svend="";
my $prevTransId="";
my $mateid;
my $transFlag=0;
my $transindex;
my $rand;
my $transFlag_doulerc;
my $process_dir;
my $fasta_len;
my @vcfArray;
my @cnvnum;
my @smalldeletionsHash;
my @deletionsHash;
my @cnvtandomHash;
my @temp;
my @temp1;
my @line;
my @maternalsnpVcf=();
my @pateralsnpVcf=();
my @snpVcf;
my @maternalcnvnum;
my @paternalcnvnum;
my @snpsHash;
my @insertionsHash;
my @inversionHash;
my @translocationsHash;
my @trans1;
my @trans2;
my @maintrans;
#my %deletionsHash;
#my %insertionsHash;
my %snpsHash;
my %snpsOutFilenames;
 my %snpsOutFilenames_copy;
#my %smallinsertionsHash;
#my %smalldeletionsHash;
#my %inversionHash;
#my %cnvtandomHash;
#my %translocationsHash;
my %translocationsTrackHash;
my %chr;
my %event_hash;
#my %chrend;

use Getopt::Long;

&Getopt::Long::GetOptions(
'refgenome_dir=s' => \$refgenome_dir,
'input_vcf_like=s' => \$input_vcf_like,
'process_dir=s' => \$process_dir,
'ins_files_dir=s' => \$ins_files_dir,
'fasta_length=i' => \$fasta_len
);
if($refgenome_dir eq "" || $ins_files_dir eq "" || $input_vcf_like eq "" || $process_dir eq "" || $fasta_len eq "")
{
	die "missing arguments\n USAGE : perl perl_genome_simulator.pl  -refgenome_dir <REF CHR DIREC> -input_vcf_like <INPUT vcf likefile> -ins_files_dir <INSERT FILES DIR> -process_dir <TEMPDIR> -fasta_length <FASTA SEQ EACH LINE LENGTH>\n";
}
$refgenome_dir =~ s/\/$//g;
$ins_files_dir =~ s/\/$//g;
$process_dir =~ s/\/$//g;
#$tempdir = "TEMPDIR";
$tempdir = $process_dir;
#creating temp directory if not exists
mkdir "$tempdir", unless -d "$tempdir";
open(LOG,">$tempdir/file.log");
die "Refgenome dir not exists\n",unless -d $refgenome_dir;
die "Insertion dir not exists\n",unless -d $ins_files_dir;
#print "\n\n\nINPUT PARAM: \n"."Reference Genome Directory\t".$refgenome_dir."\n"."Input VCF \t".$input_vcf_like."\n"."CHROMOSOME END FILE\t".$chr_end_cord."\n";
print LOG "\n\n\nINPUT PARAM: \n"."Reference Genome Directory\t".$refgenome_dir."\n"."Input VCF \t".$input_vcf_like."\nINSERTIONS FILES DIRECTORY\t$ins_files_dir"."\n PROCESSDIR\t$process_dir \n FASTA SEQ LENGTH\t$fasta_len\n";

=head
open CHREND,"$chr_end_cord" or die "no chromosome end coordinates file\n";
while(<CHREND>)
{	
	chomp($_);
	$_ =~ s/\s+/ /g;
	$_ =~ s/ /\t/g;
	@temp = split("\t",$_);
	$chrend{$temp[0]} = $temp[1];	
}
=cut

open IN_VCF,"$input_vcf_like" or die "no input like vcf file found\n";
print LOG "\n\n\nReading Input VCF like file\n";
while(<IN_VCF>)
{
	chomp($_);
	$_ =~ s/\s+/ /g;
	$_ =~ s/ /\t/g;
	#@vcfArray = split(/\t/,$_);
	#print scalar(@vcfArray);
	#print $_;	
	$linenum = $.;
	if($_ !~ m/^#/)
	{
		#print $_;
		@vcfArray = ();
		@vcfArray = split(/\t/,$_);
		#print @vcfArray;
		if(@vcfArray > 7)
		{
			$vcf_chr = $vcfArray[0];
			$vcf_pos = $vcfArray[1];
			$vcf_id = $vcfArray[2];
			$vcf_ref = $vcfArray[3];
        	$vcf_alt = $vcfArray[4];
			$vcf_qual = $vcfArray[5];
        	$vcf_filter = $vcfArray[6];
			$vcf_info = $vcfArray[7];
			#print $vcf_ref."\t".$vcf_alt."\n";
			if(length($vcf_ref) ==1 && length($vcf_alt) ==1 && $vcf_ref ne $vcf_alt)
			{
				#$snpsHash{$.."_".$vcf_chr."_".$vcf_pos} = $vcf_ref."_".$vcf_alt;
				unshift(@snpsHash,[$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,$vcf_id]);
			}
			elsif($vcf_ref !~ m/[^atgcATGC]/ && $vcf_alt !~ m/[^atgcATGC]/)
			{
				
				if(length($vcf_ref) ==1 && length($vcf_alt) > 1)
				{
					#Checking for if insertion event exists,if not then event is pushed in to array
					$temp1 = int(($vcf_pos-0.1)/$fasta_len)+1;
					$temp = "$vcf_chr"."_".$temp1;
					#print $temp."\n";
					if(exists($event_hash{$temp}))
					{
						print LOG "ignoring INSERTION event line $. as multiple events occuring on this line\n";
					}
					else
					{
						unshift(@insertionsHash,[$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,1,$vcf_id,$.]);
						#print "$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,1,$vcf_id,$.,$temp1\n";
						$event_hash{$temp} = 1;
					}
				}
				elsif(length($vcf_ref) >1 && length($vcf_alt) == 1)
                {
					#Checking for if deletion event exists,if not then event is pushed in to array
                    $temp1 = int(($vcf_pos-0.1)/$fasta_len)+1;
					$k =0;
					for($i=0;$i<1+int((length($vcf_ref)-0.1)/$fasta_len);$i++)
					{	
						$temp2 = $temp1+$i;
						$temp = "$vcf_chr"."_".$temp2;
						if(exists($event_hash{$temp}))
						{
							$k =1;
						}
					}
					if($k == 0)
					{
						unshift(@deletionsHash,[$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,1,$vcf_id,$.]);
						#print "$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,1,$vcf_id,$.,$temp1\n";
						for($i=0;$i<1+int((length($vcf_ref)-0.1)/$fasta_len);$i++)
						{		
							$temp2 = $temp1+$i;
							$temp = "$vcf_chr"."_".$temp2;
							#print $temp."\n";
							$event_hash{$temp} = 2;
							
						}
					}
					else
					{
							print LOG "ignoring DELETION event line $. as multiple events occuring on this line\n";
					}
				}
				else
				{
					 print LOG "ignoring event line $.\n";
				}
			}
			else
			{
				# variables Extraction
				$svtype = $vcf_info;
				$svtype =~s/.*SVTYPE=//;
				$svtype =~s/;.*//;
				$svlen = $vcf_info;
				$svlen =~s/.*SVLEN=//;
				$svlen =~s/;.*//;
				$svend = $vcf_info;
                $svend =~s/.*END=//;
                $svend =~s/;.*//;
				
				#checking
				if($vcf_info !~ m/SVTYPE/)
                {
                    die "Ref field or Alt field contain non nulce char and found missing SVTYPE in info field on line $.\n";
                }
				
				
				if(($vcf_alt =~ m/DEL/ || $vcf_alt =~  m/INV/) && $svend eq "")
                {
                     die "Inversion or Deletion with out END INFO field on line $.\n";
                }
				
				#adding events
				#if(length($vcf_ref) ==1 && $vcf_alt =~ m/[atgcATGC]\<*\>/ && $svtype eq "INS")
				if(length($vcf_ref) ==1 && $vcf_alt =~ m/\<*\>/ && $svtype eq "INS")
				{
					#Checking for if insertion event exists,if not then event is pushed in to array
					$ins_filename = $vcf_alt;	
					$ins_filename =~ s/.*\<//;
					$ins_filename =~ s/>.*//;
					$temp1 = int(($vcf_pos-0.1)/$fasta_len)+1;
					$temp = "$vcf_chr"."_".$temp1;
					#print $temp."\n";
					if(exists($event_hash{$temp}))
					{
						print LOG "ignoring INSERTION event line $. as multiple events occuring on this line\n";
					}
					else
					{
						unshift(@insertionsHash,[$.,$vcf_chr,$vcf_pos,$vcf_ref,$ins_filename,2,$vcf_id,$vcf_info]);
						$event_hash{$temp} = 1;
					}
				
				}
				
				elsif($vcf_alt =~ m/DEL/ && $svtype eq "DEL" && $svend > $vcf_pos)
                {
                    #Checking for if deletion event exists,if not then event is pushed in to array
                    $temp1 = int(($vcf_pos-0.1)/$fasta_len)+1;
					$k =0;
					for($i=0;$i<1+int((($svend-$vcf_pos)-0.1)/$fasta_len);$i++)
					{	
						$temp2 = $temp1+$i;
						$temp = "$vcf_chr"."_".$temp2;
						if(exists($event_hash{$temp}))
						{
							$k =1;
						}
						
					}
					if($k == 0)
					{
						unshift(@deletionsHash,[$.,$vcf_chr,$vcf_pos,$vcf_ref,$svend,2,$vcf_id,$vcf_info]);
						#print "$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,1,$vcf_id,$.,$temp1\n";
						for($i=0;$i<1+int((length($vcf_ref)-0.1)/$fasta_len);$i++)
						{		
							$temp2 = $temp1+$i;
							$temp = "$vcf_chr"."_".$temp2;
							#print $temp."\n";
							$event_hash{$temp} = 2;
						}
					}
					else
					{
							print LOG "ignoring DELETION event line $. as multiple events occuring on this line\n";
					}
				
				
				}

				elsif($vcf_alt =~ m/INV/ && $svtype eq "INV" && $svend > $vcf_pos)
				{
					#Checking for if invertion event exists,if not then event is pushed in to array
                    $temp1 = int(($vcf_pos-0.1)/$fasta_len)+1;
					$k =0;
					for($i=0;$i<1+int((($svend-$vcf_pos)-0.1)/$fasta_len);$i++)
					{	
						$temp2 = $temp1+$i;
						$temp = "$vcf_chr"."_".$temp2;
						if(exists($event_hash{$temp}))
						{
							$k =1;
						}
						
					}
					if($k == 0)
					{
						unshift(@inversionHash,[$.,$vcf_chr,$vcf_pos,$vcf_ref,$svend,$vcf_id,$vcf_info]);
						for($i=0;$i<1+int((length($vcf_ref)-0.1)/$fasta_len);$i++)
						{		
							$temp2 = $temp1+$i;
							$temp = "$vcf_chr"."_".$temp2;
							#print $temp."\n";
							$event_hash{$temp} = 3;
						}
					}
					else
					{
							print LOG "ignoring INVERTION event line $. as multiple events occuring on this line\n";
					}
				}
				
				elsif($vcf_alt =~ m/DUP:TANDEM/ && $svtype eq "DUP" && $svend > $vcf_pos)
				{
					
					#Checking for if CNV event exists,if not then event is pushed in to array
					$temp1 = int(($vcf_pos-0.1)/$fasta_len)+1;
					$k =0;
					for($i=0;$i<1+int((($svend-$vcf_pos)-0.1)/$fasta_len);$i++)
					{	
						$temp2 = $temp1+$i;
						$temp = "$vcf_chr"."_".$temp2;
						if(exists($event_hash{$temp}))
						{
							$k =1;
						}
						
					}
					if($k == 0)
					{
						unshift(@cnvtandomHash,[$.,$vcf_chr,$vcf_pos,$vcf_ref,$svend,$vcf_id,$vcf_info]);
						for($i=0;$i<1+int((length($vcf_ref)-0.1)/$fasta_len);$i++)
						{		
							$temp2 = $temp1+$i;
							$temp = "$vcf_chr"."_".$temp2;
							#print $temp."\n";
							$event_hash{$temp} = 4;
						}
					}
					else
					{
							print LOG "ignoring CNV event line $. as multiple events occuring on this line\n";
					}
				}
				elsif($svtype eq "BND" && $vcf_alt =~ m/$vcf_ref/ && ($vcf_alt =~ m/\[\d+:\d+\[/ || $vcf_alt =~ m/\]\d+:\d+\]/))
				{
					$mateid = $vcf_info;
                    $mateid =~s/.*MATEID=//;
                    $mateid =~s/;.*//;
					if($prevTransId eq "")
					{
						$prevTransId = $vcf_id;
						 @trans1 = [$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,$vcf_id,$vcf_info];
						 
					}
					else
					{
						if($prevTransId eq  $mateid)
						{
							@trans2 = [$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,$vcf_id,$vcf_info];
							push(@translocationsHash,@trans1);
							push(@translocationsHash,@trans2);
							@trans1=();
							@trans2=();
							$prevTransId = "";
						}
						else
						{
							$prevTransId = $vcf_id;
							@trans1 = [$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,$vcf_id,$vcf_info];
						}
					}
					#$translocationsHash{$.."_".$vcf_chr."_".$vcf_pos} = $vcf_ref."_".$vcf_alt;
					
					#push(@translocationsHash,@trans1 = [$.,$vcf_chr,$vcf_pos,$vcf_ref,$vcf_alt,$vcf_id,$vcf_info]);
					
				}
					
				
				else
				{
					 print LOG "ignoring event line $.\n";
				}		
				
			}
		}
		else
		{
			die "enough fields not found on line $.\n";
		}
	}
	else
	{
		print LOG "ignoring headers $.\n";
	}
}
print LOG "done reading input vcf\n\n\n";	

###########printing event hash####################
#while(($key,$value)=each %event_hash)
#{
#	print "$key\t$value\n";
#}


###########deleting event hash####################
 undef %event_hash; 
 

#print @$_, "\n" foreach ( @snpsHash );
#sorting snp events first by chr and then by coordinates
print LOG "\n\n\n\n\n\n\n#################################PROCESSING SNP EVENTS###########################\n\n";
print LOG "\nsorting SNP events\n";
@snpsHash = sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]} @snpsHash;
#print LOG "@$_\n" foreach ( @snpsHash );
############################SNPS########################
$flag =0;
print LOG "\n\ncreating SNP events for Maternal arm\n";
@maternalsnpVcf =&snp("maternal");
print LOG "@maternalsnpVcf\n";
$flag=1;
#CREATING SNPS FOR PATERNAL
print LOG "\n\ncreating SNP events for Paternal arm\n";
@paternalsnpVcf =&snp("paternal");
print LOG "@paternalsnpVcf\n";
#CREATING VCF
#$flag =0;
&createvcf;
 %snpsOutFilenames_copy = %snpsOutFilenames;

################INSERTIONS############################
@maternalsnpVcf=();
@paternalsnpVcf=();
$flag =0;
@maternalsnpVcf=&ins("maternal");
print LOG "\n\n\n\n\n\n\n#################################PROCESSING INSERTION EVENTS###########################\n\n";
print LOG "\nsorting INSERTION events\n";
#sorting small insertions events first by chr and then by coordinate
@insertionsHash = sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]} @insertionsHash;
#print LOG "@$_\n" foreach (@insertionsHash);
print LOG "\n\ncreating INSERTION events for Maternal arm\n";
@maternalsnpVcf=&ins("maternal");
print LOG "@maternalsnpVcf\n";
$flag=1;
print LOG "\n\ncreating INSERTION events for Paternal arm\n";
@paternalsnpVcf =&ins("paternal");
print LOG "@paternalsnpVcf\n";
$flag =0;

&createinsvcf;






while(my($key,$value) = each %snpsOutFilenames)
{
	if(exists($snpsOutFilenames{$key}) && exists($snpsOutFilenames_copy{$key})  && $snpsOutFilenames{$key} ne $snpsOutFilenames_copy{$key})
	{
		#print "$snpsOutFilenames{$key}\t$snpsOutFilenames_copy{$key}\n";	
		 if ((-e "$tempdir/maternal_$snpsOutFilenames_copy{$key}") && (-e "$tempdir/paternal_$snpsOutFilenames_copy{$key}")) 
		 {
			system("rm $tempdir/maternal_$snpsOutFilenames_copy{$key}\n");
			system("rm $tempdir/paternal_$snpsOutFilenames_copy{$key}\n");
		}	
		else
		{
			print LOG  "error occured in clearing files ,files not occured\n maternal_$snpsOutFilenames{$key}\n paternal_$snpsOutFilenames_copy{$key} \n";
		}
	}
}
%snpsOutFilenames_copy = %snpsOutFilenames;


#######################DELETIONS################################
print LOG "\n\n\n\n\n\n\n#################################PROCESSING DELETION EVENTS###########################\n\n";
print LOG "\nsorting DELETION events\n";
@maternalsnpVcf=();
@paternalsnpVcf=();
$flag =0;
@deletionsHash = sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]} @deletionsHash;
#print LOG "@$_\n" foreach (@deletionsHash);
print LOG "\n\ncreating DELETION events for Maternal arm\n";
@maternalsnpVcf=&del("maternal");
print LOG "@maternalsnpVcf\n";
$flag=1;
print LOG "\n\ncreating DELETION events for Paternal arm\n";
@paternalsnpVcf =&del("paternal");
print LOG "@paternalsnpVcf\n";
$flag =0;
&createdelvcf;

while(my($key,$value) = each %snpsOutFilenames)
{
	if(exists($snpsOutFilenames{$key}) && exists($snpsOutFilenames_copy{$key})  && $snpsOutFilenames{$key} ne $snpsOutFilenames_copy{$key})
	{
		#print "$snpsOutFilenames{$key}\t$snpsOutFilenames_copy{$key}\n";	
		 if ((-e "$tempdir/maternal_$snpsOutFilenames_copy{$key}") && (-e "$tempdir/paternal_$snpsOutFilenames_copy{$key}")) 
		 {
			system("rm $tempdir/maternal_$snpsOutFilenames_copy{$key}\n");
			system("rm $tempdir/paternal_$snpsOutFilenames_copy{$key}\n");
		}	
		else
		{
			print LOG  "error occured in clearing files ,files not occured\n maternal_$snpsOutFilenames{$key}\n paternal_$snpsOutFilenames_copy{$key} \n";
		}
	}
}
%snpsOutFilenames_copy = %snpsOutFilenames;


###############################CNV############################
print LOG "\n\n\n\n\n\n\n#################################PROCESSING CNV EVENTS###########################\n\n";
print LOG "\nsorting CNV events\n";
@maternalsnpVcf=();
@paternalsnpVcf=();
$flag =0;
@cnvtandomHash = sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]} @cnvtandomHash;
#print LOG "@$_\n" foreach (@cnvtandomHash);
print LOG "\n\ncreating CNV events for Maternal arm\n";
($ref1,$ref2) = cnv("maternal");
@maternalsnpVcf = @$ref1;
@maternalcnvnum = @$ref2;
print LOG "CNV events: @maternalsnpVcf\n"."CNV num: @maternalcnvnum\n";
$flag=1;
print LOG "\n\ncreating CNV events for Paternal arm\n";
($ref1,$ref2) = cnv("paternal");
@paternalsnpVcf = @$ref1;
@paternalcnvnum = @$ref2;
print LOG "CNV events: @paternalsnpVcf\n"."CNV num: @paternalcnvnum\n";
$flag =0;
&createcnvvcf;

while(my($key,$value) = each %snpsOutFilenames)
{
	if(exists($snpsOutFilenames{$key}) && exists($snpsOutFilenames_copy{$key})  && $snpsOutFilenames{$key} ne $snpsOutFilenames_copy{$key})
	{
		#print "$snpsOutFilenames{$key}\t$snpsOutFilenames_copy{$key}\n";	
		 if ((-e "$tempdir/maternal_$snpsOutFilenames_copy{$key}") && (-e "$tempdir/paternal_$snpsOutFilenames_copy{$key}")) 
		 {
			system("rm $tempdir/maternal_$snpsOutFilenames_copy{$key}\n");
			system("rm $tempdir/paternal_$snpsOutFilenames_copy{$key}\n");
		}	
		else
		{
			print LOG  "error occured in clearing files ,files not occured\n maternal_$snpsOutFilenames{$key}\n paternal_$snpsOutFilenames_copy{$key} \n";
		}
	}
}
%snpsOutFilenames_copy = %snpsOutFilenames;


###############################INVERTIONS############################
print LOG "\n\n\n\n\n\n\n#################################PROCESSING INVERTIONS EVENTS###########################\n\n";
print LOG "\nsorting INSERTIONS events\n";
@maternalsnpVcf=();
@paternalsnpVcf=();
$flag =0;
@inversionHash = sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]} @inversionHash;
#print LOG "@$_\n" foreach (@inversionHash);
print LOG "\n\ncreating INVERTION events for Maternal arm\n";
@maternalsnpVcf = inv("maternal");
print LOG "@maternalsnpVcf\n";
$flag=1;
print LOG "\n\ncreating INVERTION events for Paternal arm\n";
@paternalsnpVcf = inv("paternal");
print LOG "@paternalsnpVcf\n";
$flag =0;
&createinvvcf;

while(my($key,$value) = each %snpsOutFilenames)
{
	if(exists($snpsOutFilenames{$key}) && exists($snpsOutFilenames_copy{$key})  && $snpsOutFilenames{$key} ne $snpsOutFilenames_copy{$key})
	{
		#print "$snpsOutFilenames{$key}\t$snpsOutFilenames_copy{$key}\n";	
		 if ((-e "$tempdir/maternal_$snpsOutFilenames_copy{$key}") && (-e "$tempdir/paternal_$snpsOutFilenames_copy{$key}")) 
		 {
			system("rm $tempdir/maternal_$snpsOutFilenames_copy{$key}\n");
			system("rm $tempdir/paternal_$snpsOutFilenames_copy{$key}\n");
		}	
		else
		{
			print LOG  "error occured in clearing files ,files not occured\n maternal_$snpsOutFilenames{$key}\n paternal_$snpsOutFilenames_copy{$key} \n";
		}
	}
}
%snpsOutFilenames_copy = %snpsOutFilenames;

##############################TRANSLOCATIONS############################
print LOG "\n\n\n\n\n\n\n#################################PROCESSING TRANSLOCATIONS EVENTS###########################\n\n";
print LOG "\nTRANSLOCATIONS events\n";
@maternalsnpVcf=();
@paternalsnpVcf=();
$flag =0;
#@translocationsHash = sort {$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2]} @translocationsHash;
#print LOG "@$_\n" foreach (@translocationsHash);
print LOG "\n\ncreating TRANSLOCATION events for Maternal arm\n";
@maternalsnpVcf = trans("maternal");
print LOG "@maternalsnpVcf\n";

#RESET HASH $translocationsTrackHash ############
#while(($key,$value)=each %translocationsTrackHash)
#{
#	delete($translocationsTrackHash{$key});
#}
foreach $key (keys %translocationsTrackHash) 
{
	delete $translocationsTrackHash{$key};
}
$flag=1;
print LOG "\n\ncreating TRANSLOCATION events for Paternal arm\n";
@paternalsnpVcf = trans("paternal");
print LOG "@paternalsnpVcf\n";
$flag =0;
&createtransvcf;


while(my($key,$value) = each %translocationsTrackHash)
{
	if(exists($translocationsTrackHash{$key}) && exists($snpsOutFilenames_copy{$key}))
	{
		#print "$translocationsTrackHash{$key}\t$snpsOutFilenames_copy{$key}\n";	
		 if ((-e "$tempdir/maternal_$snpsOutFilenames_copy{$key}") && (-e "$tempdir/paternal_$snpsOutFilenames_copy{$key}") ) 
		 {
			system("rm $tempdir/maternal_$snpsOutFilenames_copy{$key}\n");
			system("rm $tempdir/paternal_$snpsOutFilenames_copy{$key}\n");
		}	
		else
		{
			print LOG  "error occured in clearing files ,files not occured\n maternal_$snpsOutFilenames{$key}\n paternal_$snpsOutFilenames_copy{$key} \n";
		}
	}
}

opendir(DIR, $tempdir) || die "no directory exists $tempdir$!";
while($line = readdir(DIR))
{
		if($line =~ m/\.fa$/)
		{
			#print $line."\n";
			open BUFF,"$tempdir/$line" or die "no file exists $tempdir/$line\n";
			open WRBUFF,">$tempdir/tempfile" or die "unable to write to $tempdir/tempfile\n";
			$temp = <BUFF>;
			#print "$line\n";
			print WRBUFF $temp;
			$temp1 = 0;
			while($temp=<BUFF>)
			{
				chomp($temp);
				$temp =~ s/ //g;
				if($temp ne "")
				{
					print WRBUFF "$temp\n";
				}
			#	@array = ();
			#	@array = split('',$temp);
			#	for($i=0;$i<@array;$i++)
			#	{
			#		if($temp1 == $fasta_len)
			#		{
			#			print WRBUFF "\n";
			#			$temp1 = 0;
			#		}
			#		print WRBUFF $array[$i];
			#		$temp1++;
			#	}
				
			}
			#print WRBUFF "\n";
			#die;
			system("mv $tempdir/tempfile $tempdir/$line");
			
		}
}

print "Process Completed\n";

##################SUBROUTINES######################
sub snp {
 my($val) = shift @_;
$temp = "";
@snpVcf=();
foreach (@snpsHash )
{
	#print $$_[1]."\n";
#for every new chr creating file
	if($temp ne $$_[1])
	{
		if($temp ne "")
		{

			print $maternal $line."\n";
			#print $paternal $line."\n";
			while($line=<$buff>)
			{
				print $maternal $line;
				#print $paternal $line;
			}
		}
		$buff = "CHR$$_[1]";
		$maternal = "MATERNAL$$_[1]";
		#$paternal = "PATERNAL$$_[1]";
		#$file = "$refgenome_dir/$$_[1]".'.fa';
		$file = "$refgenome_dir/chr$$_[1]".'.fa';
		open $buff,"$file" or die "$file\n";	
		$file = "SNPchr$$_[1]".'.fa';
		$snpsOutFilenames{$$_[1]} =  $file;
		$file = "$val"."_SNPchr$$_[1]".'.fa';
		
		#opening the ref file
		open($maternal,">$tempdir/$file");
		#$file = "PATERNAL_SNPchr$$_[1]".'.fa';
		#open($paternal,">$tempdir/$file");
		$line=<$buff>;
		print $maternal $line;
		#print $paternal $line;	
		$length =0;
		$prevlength =0;
		#print $line."\n";
	}
	
	#extracting the line and checking if coordinate mentioned more than size of chr
	while($$_[2] > $length)
	{
		if($length != 0)
		{
			print $maternal $line."\n";
			#print $paternal $line;
		}
		#print "sucess\n";
		if($line=<$buff>)
		{
	
			chomp($line);
			$prevlength = $length;
			$length=$length+length($line);
		}
		else
		{
			#push(@snpVcf,2);
			die "event $$_[0] coord exceeds chr length please change input vcf file and rerun ignoring event not possible (End of chr file)\n";
			
		}
	} 	
	#print $$_[2]."\t".$length."\n";
	#die;
	@temp=();
	@temp = split('',$line);
	#print $line."\n";
	#$v = $$_[2]-$prevlength-1;
		
	#print $v."\t".$prevlength."\t".$line."\t".$temp[$$_[2]-$prevlength-1]."\t".$$_[3]."\t".$$_[1]."\n";

		#checking if reference if equal is equal
	if(uc($temp[$$_[2]-$prevlength-1]) eq uc($$_[3]))
	{
			#print "testnum1 $maternalsnpVcf[@snpVcf]\n";
			if(rand() < 0.5)
			{	
				$temp[$$_[2]-$prevlength-1] = $$_[4];
				$line=join('',@temp);
				push(@snpVcf,1);
				#print @snpVcf."\n";
			}	
			else
			{
				 #print @snpVcf."\n";
				#print $flag."\t".$maternalsnpVcf[@paternalsnpVcf]."\n";
				if($flag ==1 && $maternalsnpVcf[@snpVcf] == 0)
				{
					$temp[$$_[2]-$prevlength-1] = $$_[4];
					$line=join('',@temp);
					push(@snpVcf,1);
				}  
				else
				{
					$line=join('',@temp);
					push(@snpVcf,0);
				}
			}
		
	}
	else
	{
		print LOG "ignoring event line $$_[0] snp reference not matched \n";
		#print LOG "\t ".uc($temp[$$_[2]-$prevlength-1])." eq ".uc($$_[3])."\n";
		push(@snpVcf,2);
	}
	$temp = $$_[1];
	
}
if($temp ne "")
{
	print $maternal $line."\n";
        #print $paternal $line."\n";
        while($line=<$buff>)
        {
        	print $maternal $line;
            #print $paternal $line;
        }
}

close(my $buff);
close($maternal);
return(@snpVcf);
}

#####################subroutine create vcf###################
sub createvcf{
$index = 0;

open(VCF,">$tempdir/file.vcf");
print VCF '##fileformat=VCFv4.1';
print VCF "\n";
print VCF '##fileDate=20090805';
print VCF "\n";
print VCF '##source=myImputationProgramV3.1';
print VCF "\n";
print VCF '##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta';
print VCF "\n";
print VCF '##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>';
print VCF "\n";
print VCF '##phasing=partial';
print VCF "\n";
print VCF '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">';
print VCF "\n";
print VCF '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">';
print VCF "\n";
print VCF '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">';
print VCF "\n";
print VCF '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">';
print VCF "\n";
print VCF '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">';
print VCF "\n";
print VCF '##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">';
print VCF "\n";
print VCF '##FILTER=<ID=q10,Description="Quality below 10">';
print VCF "\n";
print VCF '##FILTER=<ID=s50,Description="Less than 50% of samples have data">';
print VCF "\n";
print VCF '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">';
print VCF "\n";
print VCF '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">';
print VCF "\n";
print VCF '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">';
print VCF "\n";
print VCF '##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">';
print VCF "\n";
print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE";
print VCF "\n";	
	
	for(@snpsHash )
	{
		#print $index."\n";
		if($maternalsnpVcf[$index] <2 && $paternalsnpVcf[$index] <2)
		{
			#if($paternalsnpVcf[$index] < $maternalsnpVcf[$index])
			#{
			#	$val = $paternalsnpVcf[$index]."\/".$maternalsnpVcf[$index];
			#}
			#else
			#{
			#	$val = $maternalsnpVcf[$index]."\/".$paternalsnpVcf[$index];
			#}
			$val = $maternalsnpVcf[$index]."\|".$paternalsnpVcf[$index];
			$file = "$$_[1]\t$$_[2]\t$$_[5]\t$$_[3]\t$$_[4]\t"."\.\t"."PASS\t\.\tGT\t$val\n"; 
			print VCF $file;
		}
		$index++;
	}
	close(VCF);
}

########################################CREATE VCF INSERTIONS########################
sub createinsvcf{
	#my($val) = shift @_;
	#print "$val\n@{$val}\n";
	$index = 0;
	open(VCF,">>$tempdir/file.vcf");
	for(@insertionsHash )
	#for(@{"$val"})
        {
                #print $index."\n";
                if($maternalsnpVcf[$index] <2 && $paternalsnpVcf[$index] <2)
                {
                        #if($paternalsnpVcf[$index] < $maternalsnpVcf[$index])
                        #{
                        #        $val = $paternalsnpVcf[$index]."\/".$maternalsnpVcf[$index];
                        #}
                        #else
                        #{
                        #        $val = $maternalsnpVcf[$index]."\/".$paternalsnpVcf[$index];
                        #}
						$val = $maternalsnpVcf[$index]."\|".$paternalsnpVcf[$index];
						if($$_[5] == 2)
						{
							$file = "$$_[1]\t$$_[2]\t$$_[6]\t$$_[3]\t".'<'.$$_[4].'>'."\t"."\.\t"."PASS\t$$_[7]\tGT\t$val\n";
                        }
						else
						{
							$file = "$$_[1]\t$$_[2]\t$$_[6]\t$$_[3]\t$$_[4]\t"."\.\t"."PASS\t\.\tGT\t$val\n";
						}
						print VCF $file;
                }
                $index++;
        }
        close(VCF);
}

########################################CREATE VCF DELETIONS########################
sub createdelvcf{
	#my($val) = shift @_;
	#print "$val\n@{$val}\n";
	$index = 0;
	open(VCF,">>$tempdir/file.vcf");
	for(@deletionsHash )
	#for(@{"$val"})
        {
                #print $index."\n";
                if($maternalsnpVcf[$index] <2 && $paternalsnpVcf[$index] <2)
                {
                       # if($paternalsnpVcf[$index] < $maternalsnpVcf[$index])
                        #{
                        #        $val = $paternalsnpVcf[$index]."\/".$maternalsnpVcf[$index];
                        #}
                        #else
                        #{
                        #        $val = $maternalsnpVcf[$index]."\/".$paternalsnpVcf[$index];
                        #}
						$val = $maternalsnpVcf[$index]."\|".$paternalsnpVcf[$index];
						if($$_[5] == 2)
						{	
							if($maternalsnpVcf[$index] == 0 || $paternalsnpVcf[$index] == 0)
							{
								$file = "$$_[1]\t$$_[2]\t$$_[6]\t$$_[3]\t<DEL>\t"."\.\t"."PASS\t$$_[7]\tGT:CN\t$val".":1\n";
							}
							else
							{
								$file = "$$_[1]\t$$_[2]\t$$_[6]\t$$_[3]\t<DEL>\t"."\.\t"."PASS\t$$_[7]\tGT:CN\t$val".":0\n";
							}
						}
						else
						{
							$file = "$$_[1]\t$$_[2]\t$$_[6]\t$$_[3]\t$$_[4]\t"."\.\t"."PASS\t\.\tGT\t$val\n";
						}
						print VCF $file;
                }
                $index++;
        }
        close(VCF);
}

########################################CREATE VCF INVERTIONS########################
sub createinvvcf{
	#my($val) = shift @_;
	#print "$val\n@{$val}\n";
	$index = 0;
	open(VCF,">>$tempdir/file.vcf");
	for(@inversionHash)
	{
                #print $index."\n";
                if($maternalsnpVcf[$index] <2 && $paternalsnpVcf[$index] <2)
                {
                       # if($paternalsnpVcf[$index] < $maternalsnpVcf[$index])
                       # {
                       #         $val = $paternalsnpVcf[$index]."\/".$maternalsnpVcf[$index];
                       # }
                       # else
                       # {
                       #         $val = $maternalsnpVcf[$index]."\/".$paternalsnpVcf[$index];
                       # }
						$val = $maternalsnpVcf[$index]."\|".$paternalsnpVcf[$index];
						$file = "$$_[1]\t$$_[2]\t$$_[5]\t$$_[3]\t<INV>\t"."\.\t"."PASS\t$$_[6]\tGT\t$val\n";
						
						#print $file."test\n";
						print VCF $file;
                }
                $index++;
        }
        close(VCF);
}

########################################CREATE VCF CNV########################
sub createcnvvcf{
	#my($val) = shift @_;
	#print "$val\n@{$val}\n";
	$index = 0;
	open(VCF,">>$tempdir/file.vcf");
	for(@cnvtandomHash )
	#for(@{"$val"})
        {
                #print $index."\n";
                if($maternalsnpVcf[$index] <2 && $paternalsnpVcf[$index] <2)
                {
                       # if($paternalsnpVcf[$index] < $maternalsnpVcf[$index])
                       # {
                       #         $val = $paternalsnpVcf[$index]."\/".$maternalsnpVcf[$index];
                       # }
                       # else
                       # {
                       #         $val = $maternalsnpVcf[$index]."\/".$paternalsnpVcf[$index];
                       # }
						$val = $maternalsnpVcf[$index]."\|".$paternalsnpVcf[$index];
						$file = "$$_[1]\t$$_[2]\t$$_[5]\t$$_[3]\t<DUP:TANDEM>\t"."\.\t"."PASS\t$$_[6]\tGT:CN\t$val";
						$val = $paternalcnvnum[$index]+$maternalcnvnum[$index]+2;
						$file = $file.":$val";
						#print $file."test\n";
						print VCF $file."\n";
                }
                $index++;
        }
        close(VCF);
}

########################################CREATE VCF TRANSLOCATION########################
sub createtransvcf{
	$transindex = -1;
	$index = 0;
	open(VCF,">>$tempdir/file.vcf");
	for($index=0;$index<@translocationsHash;$index++)
	{		
				$transindex++;	
                if($maternalsnpVcf[$transindex] <2 && $paternalsnpVcf[$transindex] <2)
                {
                        #if($paternalsnpVcf[$transindex] < $maternalsnpVcf[$transindex])
                        #{
                        #        $val = $paternalsnpVcf[$transindex]."\/".$maternalsnpVcf[$transindex];
                        #}
                        #else
                        #{
                        #        $val = $maternalsnpVcf[$transindex]."\/".$paternalsnpVcf[$transindex];
                        #}
						#$file = ${$translocationsHash[$transindex]}[1];
						#print $file."\n";
						$val = $maternalsnpVcf[$transindex]."\|".$paternalsnpVcf[$transindex];
						$file = "${$translocationsHash[$index]}[1]\t${$translocationsHash[$index]}[2]\t${$translocationsHash[$index]}[5]\t${$translocationsHash[$index]}[3]\t${$translocationsHash[$index]}[4]\t"."\.\t"."PASS\t${$translocationsHash[$index]}[6]\tGT\t$val";
						#print $file."\n";
						print VCF $file."\n";
						$index++;
						$file = "${$translocationsHash[$index]}[1]\t${$translocationsHash[$index]}[2]\t${$translocationsHash[$index]}[5]\t${$translocationsHash[$index]}[3]\t${$translocationsHash[$index]}[4]\t"."\.\t"."PASS\t${$translocationsHash[$index]}[6]\tGT\t$val";
						#print $file."\n";
						print VCF $file."\n";
						
                }
                else
				{
					$index++;
				}
	}
        close(VCF);
}

#######################################INSERTION SUBROUTINE###########################
sub ins{
$temp = "";
my($val) = shift @_;
@snpVcf=();
foreach (@insertionsHash)
{
	if($temp ne $$_[1])
	{
		$v1 = 0;
		if($temp ne "")
		{

			print $maternal $line."\n";
			while($line=<$buff>)
			{
				print $maternal $line;
			}
		}
		#buffers
		$buff = "CHR$$_[1]";
		$maternal = "MATERNAL$$_[1]";
		if(exists($snpsOutFilenames{$$_[1]}))
		{
			$file = "$tempdir/$val"."_SNPchr$$_[1]".'.fa';
		}
		else
		{
			$file = "$refgenome_dir/chr$$_[1]".'.fa';
		}
		#print $file."test\n";
		open $buff,"$file" or die "$file\n";
		
		$file = "INSchr$$_[1]".'.fa';
		if($flag== 1)
		{
			$snpsOutFilenames{$$_[1]} =  $file;
		}
		$file = "$val"."_INSchr$$_[1]".'.fa';
		open($maternal,">$tempdir/$file");
		$line=<$buff>;
		print $maternal $line;
		$length =0;
		$prevlength =0;
	}
	while($$_[2] > $length)
	{
		$v1 =0;
		if($length != 0)
		{
			print $maternal $line."\n";
			
		}
		
		if($line=<$buff>)
		{
	
			chomp($line);
			$prevlength = $length;
			$length=$length+length($line);
		}
		else
		{
			
			die "event $$_[0] coord exceeds chr length please change input vcf file and rerun ignoring event not possible (End of chr file)\n";
			
		}
	} 	
	#print $$_[2]."\t".$length."\n";
	#die;
	@temp=();
	@temp = split('',$line);
	#print "$$_[0]\t".scalar(@temp)."\n";
	#print $line."\n";
	#$v = $$_[2]-$prevlength-1;
	#print $v."\t".$prevlength."\t".$line."\n";	
	#print $v."\t".$prevlength."\t".$line."\t".$temp[$$_[2]-$prevlength-1]."\t".$$_[3]."\t".$$_[1]."\n";	
	#if(uc($temp[$$_[2]-$prevlength-1]) eq uc($$_[3]))
	#{
			if($$_[5] == 2)
			{
				open INSFILE,"$ins_files_dir/$$_[4]" or die "no file exists $$_[4] ";
				$line=<INSFILE>;
				$line="";
				@line = <INSFILE>;
				#$$_[4]=join('',@line);
				$temp2 = join('',@line);
				@line=();
				#$$_[4] =~ s/\n//g;
				#$$_[4] =$$_[3].$$_[4];
				$temp2 =~ s/\n//g;
				$temp2 =$$_[3].$temp2;
				close(INSFILE);
				
			}
			else
			{

				$temp2 = $$_[4];
			}
			#print $temp2."\ttest\n";
			#$v1 = @temp -($length-$prevlength);
		
			#print $v1."\t".@temp."\t".$length."\t".$prevlength."\n";
			#print "testnum1 $maternalsnpVcf[@snpVcf]\n";
			#print "maintest\t".$temp[$$_[2]-$prevlength-1]."\t".$$_[3]."\n";
			#die;
			#$exam = $$_[2]-$prevlength-1;
			#print "maintest $exam @temp\n";
			#if(uc($temp[$$_[2]-$prevlength-1]) ne uc($$_[3]) && $v1 == 0)
			#{
			#	$v1 =2;
			#}
			if(rand() < 0.5 && $v1 ==0)
			{
				$temp[$$_[2]-$prevlength-1] = $temp2;
				$line=join('',@temp);
				push(@snpVcf,1);
				#print @snpVcf."\n";
			}	
			else
			{
				 #print @snpVcf."\n";
				#print "$v1\t$flag"."\t".$maternalsnpVcf[@snpVcf]."\n";
				#if($v1 == 2)
				#{
				#	print "ignore event $$_[0] as reference base $$_[3] not matched file ref $temp[$$_[2]-$prevlength-1] \n";
				#	$line=join('',@temp);
				#	push(@snpVcf,2);
				#}
				#elsif($v1 ==1)
				if($v1 !=0)
				{
					#print LOG "ignore event $$_[0] multiple events on the same line\n";
					$line=join('',@temp);
					push(@snpVcf,2);
				}
				elsif($flag ==1 && $maternalsnpVcf[@snpVcf] == 0)
				{
				#	print "sucess\n";
					$temp[$$_[2]-$prevlength-1] = $temp2;
					$line=join('',@temp);
					push(@snpVcf,1);
				}  
				else
				{
					$line=join('',@temp);
					push(@snpVcf,0);
				}
			}
	#}
	#else
	#{
	#	print "ignoring event line $$_[0] snp reference not matched \n";
	#	push(@snpVcf,2);
	#}
	$v1 =1;
	$temp = $$_[1];		
}
if($temp ne "")
{
	print $maternal $line."\n";
        #print $paternal $line."\n";
        while($line=<$buff>)
        {
        	print $maternal $line;
            #print $paternal $line;
        }
}
close($buff);
close($maternal);
return(@snpVcf);
}

############################################DELETION##########################################
sub del{
$temp = "";
$index = 1;
my($val) = shift @_;
@snpVcf=();
foreach (@deletionsHash)
{
	$v1 = 0;
	if($temp ne $$_[1])
	{
		
		if($temp ne "")
		{

			print $maternal $line."\n";
			while($line=<$buff>)
			{
				print $maternal $line;
			}
		}
		#buffers
		$buff = "CHR$$_[1]";
		$maternal = "MATERNAL$$_[1]";
		if(exists($snpsOutFilenames{$$_[1]}))
		{
			
			$file = "$tempdir/$val"."_$snpsOutFilenames{$$_[1]}";
		}
		else
		{
			$file = "$refgenome_dir/chr$$_[1]".'.fa';
		}
		
		open $buff,"$file" or die "$file\n";
		
		$file = "DELchr$$_[1]".'.fa';
		if($flag== 1)
		{
			$snpsOutFilenames{$$_[1]} =  $file;
		}
		$file = "$val"."_DELchr$$_[1]".'.fa';
		open($maternal,">$tempdir/$file");
		$line=<$buff>;
		print $maternal $line;
		$length =0;
		$prevlength =0;
	}
	if($$_[2] <= $length)
	{
		$v1 =1;
	}
	while($$_[2] > $length)
	{
		
		if($length != 0)
		{
			print $maternal $line."\n";
			
		}
		
		if($line=<$buff>)
		{
	
			chomp($line);
			$prevlength = $length;
			#$length=$length+length($line);
			$length=$length+$fasta_len;
		}
		else
		{
			
			die "event $$_[0] coord exceeds chr length please change input vcf file and rerun ignoring event not possible (End of chr file)\n";
			
		}
	} 	
		
	


			$eventstart =  $prevlength;
			if($$_[5] == 1)
			{
				
				$eventstop =length($$_[3]) -1+$$_[2];
				
				
			}
			else
			{

				$eventstop = $$_[4];
			}
			$counter = 1;
			$temp2 = $line;
			$temp_fail = $line."\n";
			while($length  < $eventstop)
			{	$line=<$buff>;
				chomp($line);
				$prevlength = $length;
				$counter++;
				$length=$length+$fasta_len;
				$temp2 = $temp2.$line;	
				$temp_fail = $temp_fail.$line."\n";
			}
			if(length($temp2) != $counter*$fasta_len)
			{
				$v1 = 1;
			}
			$eventstop = $eventstop-$eventstart-1;
			$eventstart = $$_[2]-$eventstart-1;
			
			#print "counter $counter\t prevlength $prevlength\tlength $length\ttemp ".length($temp2)."\tv1 $v1\t eventstart $eventstart\t eventstop $eventstop\n"; 
			
			#print "$temp2\n";
			

			#print "$$_[3]\t$$_[2]\t".$temp2."\ttest\n";
			#$v1 = @temp -($length-$prevlength);
		
			#print $v1."\t".@temp."\t".$length."\t".$prevlength."\n";
			#print "testnum1 $maternalsnpVcf[@snpVcf]\n";
			if($flag ==1 && $maternalsnpVcf[@snpVcf] == 2)
			{
				 $v1  = 1;
			}
			if(rand() < 0.5 && $v1 ==0)
			{
				#print "test1 $temp2\n";
				@temp = split('',$temp2);
				for($i=0;$i<@temp;$i++)
				{
					if($i>$eventstart && $i<=$eventstop)
					{
						$temp[$i] = " ";
					}
					if((($i+1)%$fasta_len) == 0)
					{
						#print "success\n";
						$temp[$i] = $temp[$i]."\n";
					}
				}
				$temp2 = join('',@temp);
				$temp2 =~ s/\n$//g;
				#print "test2 $temp2\n";
				$line=$temp2;
				push(@snpVcf,1);
				#print @snpVcf."\n";
			}	
			else
			{
				 #print @snpVcf."\n";
				#print "$v1\t$flag"."\t".$maternalsnpVcf[@snpVcf]."\n";
				if($v1 !=0)
				{
					print "ignore event $$_[0] multiple events occuring on the same lines of $fasta_len bp length\n";
					$temp_fail =~ s/\n$//g;
					$line=$temp_fail;	
					push(@snpVcf,2);
				}
				elsif($flag ==1 && $maternalsnpVcf[@snpVcf] == 0)
				{
					#print "sucess $$_[0] $maternalsnpVcf[@snpVcf]\n";
					@temp = split('',$temp2);
					for($i=0;$i<@temp;$i++)
					{
						if($i>$eventstart && $i<=$eventstop)
						{
							$temp[$i] = " ";
						}
						if((($i+1)%$fasta_len) == 0)
						{
							#print "success\n";
							$temp[$i] = $temp[$i]."\n";
						}
					}
					$temp2 = join('',@temp);
					$temp2 =~ s/\n$//g;
					$line=$temp2;
					push(@snpVcf,1);
				}  
				else
				{	
					$temp_fail =~ s/\n$//g;
					$line=$temp_fail;
					push(@snpVcf,0);
				}
			}

	
	$temp = $$_[1];		
	$index++;
}

if($temp ne "")
{
	print $maternal $line."\n";
        #print $paternal $line."\n";
        while($line=<$buff>)
        {
        	print $maternal $line;
            #print $paternal $line;
        }
}

close($buff);
close($maternal);
return(@snpVcf);

}

############################################CNV##########################################
sub cnv{
$temp = "";
$index = 1;
my($val) = shift @_;
@snpVcf=();
@cnvnum=();
foreach (@cnvtandomHash)
{
	$v1 = 0;
	if($temp ne $$_[1])
	{
		
		if($temp ne "")
		{

			print $maternal $line."\n";
			while($line=<$buff>)
			{
				print $maternal $line;
			}
		}
		#buffers
		$buff = "CHR$$_[1]";
		$maternal = "MATERNAL$$_[1]";
		if(exists($snpsOutFilenames{$$_[1]}))
		{
			
			$file = "$tempdir/$val"."_$snpsOutFilenames{$$_[1]}";
		}
		else
		{
			$file = "$refgenome_dir/chr$$_[1]".'.fa';
		}
		#print "test $buff $maternal $file\n";
		open $buff,"$file" or die " no file exists $file\n";
		
		$file = "CNVchr$$_[1]".'.fa';
		if($flag== 1)
		{
			$snpsOutFilenames{$$_[1]} =  $file;
		}
		$file = "$val"."_CNVchr$$_[1]".'.fa';
		open($maternal,">$tempdir/$file");
		$line=<$buff>;
		#print "test 2 $line\n";
		print $maternal $line;
		$length =0;
		$prevlength =0;
	}

	if($$_[2] <= $length)
	{
		$v1 =1;
	}
	while($$_[2] > $length)
	{
		
		if($length != 0)
		{
			print $maternal $line."\n";
			
		}
		
		if($line=<$buff>)
		{
	
			chomp($line);
			$prevlength = $length;
			#$length=$length+length($line);
			$length=$length+$fasta_len;
		}
		else
		{
			
			die "event $$_[0] coord exceeds chr length please change input vcf file and rerun ignoring event not possible (End of chr file)\n";
			
		}
	} 	
	


			$eventstart =  $prevlength;
			$eventstop = $$_[4];
			$counter = 1;
			$temp2 = $line;
			$temp_fail = $line."\n";
			while($length  < $eventstop)
			{	$line=<$buff>;
				chomp($line);
				$prevlength = $length;
				$counter++;
				$length=$length+$fasta_len;
				$temp2 = $temp2.$line;
				$temp_fail = $temp_fail."$line\n";
			}
			$temp3 = $temp2;
			$temp3 =~ s/ //g;
			if(length($temp3) != $counter*$fasta_len)
			{
				$v1 = 1;
			}
			$eventstop = $eventstop-$eventstart-1;
			$eventstart = $$_[2]-$eventstart-1;
			
			#print "counter $counter\t prevlength $prevlength\tlength $length\ttemp ".length($temp2)."\tv1 $v1\t eventstart $eventstart\t eventstop $eventstop\n"; 
			
		#	print "$temp2\n";
			
		#	print $maincnv."\t$maincnvnum\n$temp2\n";
			


			#print "$$_[3]\t$$_[2]\t".$temp2."\ttest\n";
			#$v1 = @temp -($length-$prevlength);
		
			#print $v1."\t".@temp."\t".$length."\t".$prevlength."\n";
			#print "testnum1 $maternalsnpVcf[@snpVcf]\n";
			if($flag ==1 && $maternalsnpVcf[@snpVcf] == 2)
			{
				 $v1  = 1;
			}
			$flag2 = 0;
			if($flag ==1 && $maternalsnpVcf[@snpVcf] == 1)
			{
				$flag2 = 1;
			}	
			if(rand() < 0.5 && $v1 ==0 && $flag2 == 0)
			{
				#print "test1 $temp2\n";
				@temp = split('',$temp2);
				$cnv="";
				for($i=0;$i<@temp;$i++)
				{
					if($i>$eventstart && $i<=$eventstop)
					{
						$cnv = $cnv.$temp[$i];
					}
					if((($i+1)%$fasta_len) == 0)
					{
						$temp[$i] = $temp[$i]."\n";
					}
				}
			
				$temp2 = join('',@temp);
				$temp2 =~ s/\n$//g;
				$maincnvnum = int(rand(9)-0.1)+2;
				$maincnv = "";
				for($i=0;$i<$maincnvnum;$i++)
				{
					$maincnv=$maincnv.$cnv;
				}
				@temp = split('',$temp2);
				$temp[$eventstart] = $temp[$eventstart].$maincnv;
				$temp2 = join('',@temp);
				#$temp2 = $maincnv.$temp2;
				#print "$maincnvnum\t$maincnv\n";
				$line=$temp2;
				push(@snpVcf,1);
				push(@cnvnum,$maincnvnum);
				#print @snpVcf."\n";
			}	
			else
			{
				 #print @snpVcf."\n";
				#print "$v1\t$flag"."\t".$maternalsnpVcf[@snpVcf]."\n";
				if($v1 !=0)
				{
					print LOG "ignore event $$_[0] multiple events occuring on the same lines of $fasta_len bp length\n";
					$temp_fail =~ s/\n$//g;
					$line=$temp_fail;	
					push(@snpVcf,2);
					push(@cnvnum,0);
				}

				elsif($flag ==1 && $maternalsnpVcf[@snpVcf] == 0)
				{
					#print "sucess $$_[0] $maternalsnpVcf[@snpVcf]\n";
					@temp = split('',$temp2);
					$cnv="";
					for($i=0;$i<@temp;$i++)
					{
						if($i>$eventstart && $i<=$eventstop)
						{
							$cnv = $cnv.$temp[$i];
						}
						if((($i+1)%$fasta_len) == 0)
						{
							$temp[$i] = $temp[$i]."\n";
						}
					}
			
					$temp2 = join('',@temp);
					$temp2 =~ s/\n$//g;
					$maincnvnum = int((rand(9)-0.1))+2;
					$maincnv = "";
					for($i=0;$i<$maincnvnum;$i++)
					{
						$maincnv=$maincnv.$cnv;
					}
					$temp[$eventstart] = $temp[$eventstart].$maincnv;
					$temp2 = join('',@temp);
					#$temp2 = $maincnv.$temp2;
					#print "test2 $temp2\n";
					$line=$temp2;
					push(@snpVcf,1);
					push(@cnvnum,$maincnvnum);
				}  

				else
				{	
					$temp_fail =~ s/\n$//g;
					$line=$temp_fail;
					push(@snpVcf,0);
					push(@cnvnum,0);
				}
			}

	
	$temp = $$_[1];		
	$index++;
}

if($temp ne "")
{
	print $maternal $line."\n";
        #print $paternal $line."\n";
        while($line=<$buff>)
        {
        	print $maternal $line;
            #print $paternal $line;
        }
}

close($buff);
close($maternal);
return(\@snpVcf,\@cnvnum);

}

############################################INVERSION##########################################
sub inv{
$temp = "";
$index = 1;
my($val) = shift @_;
@snpVcf=();
foreach (@inversionHash)
{
	$v1 = 0;
	if($temp ne $$_[1])
	{
		
		if($temp ne "")
		{

			print $maternal $line."\n";
			while($line=<$buff>)
			{
				print $maternal $line;
			}
		}
		#buffers
		$buff = "CHR$$_[1]";
		$maternal = "MATERNAL$$_[1]";
		if(exists($snpsOutFilenames{$$_[1]}))
		{
			
			$file = "$tempdir/$val"."_$snpsOutFilenames{$$_[1]}";
		}
		else
		{
			$file = "$refgenome_dir/chr$$_[1]".'.fa';
		}
		
		open $buff,"$file" or die "$file\n";
		
		$file = "INVchr$$_[1]".'.fa';
		if($flag== 1)
		{
			$snpsOutFilenames{$$_[1]} =  $file;
		}
		$file = "$val"."_INVchr$$_[1]".'.fa';
		open($maternal,">$tempdir/$file");
		$line=<$buff>;
		print $maternal $line;
		$length =0;
		$prevlength =0;
	}
	if($$_[2] <= $length)
	{
		$v1 =1;
	}
	while($$_[2] > $length)
	{
		
		if($length != 0)
		{
			print $maternal $line."\n";
			
		}
		
		if($line=<$buff>)
		{
	
			chomp($line);
			$prevlength = $length;
			#$length=$length+length($line);
			$length=$length+$fasta_len;
		}
		else
		{
			
			die "event $$_[0] coord exceeds chr length please change input vcf file and rerun ignoring event not possible (End of chr file)\n";
			
		}
	} 	
		
	


			$eventstart =  $prevlength;
			$eventstop = $$_[4];
			
			$counter = 1;
			$temp2 = $line;
			$temp_fail = $line."\n";
			while($length  < $eventstop)
			{	$line=<$buff>;
				chomp($line);
				$prevlength = $length;
				$counter++;
				$length=$length+$fasta_len;
				$temp2 = $temp2.$line;	
				$temp_fail = $temp_fail.$line."\n";
			}
			if(length($temp2) != $counter*$fasta_len)
			{
				$v1 = 1;
			}
			$eventstop = $eventstop-$eventstart-1;
			$eventstart = $$_[2]-$eventstart-1;

				#print "counter $counter\t prevlength $prevlength\tlength $length\ttemp ".length($temp2)."\tv1 $v1\t eventstart $eventstart\t eventstop $eventstop\n"; 
			
				#print "entering $temp2\n";
				
				#print $temp2."\n";
		

			#print "$$_[3]\t$$_[2]\t".$temp2."\ttest\n";
			#$v1 = @temp -($length-$prevlength);
		
			#print $v1."\t".@temp."\t".$length."\t".$prevlength."\n";
			#print "testnum1 $maternalsnpVcf[@snpVcf]\n";
			if($flag ==1 && $maternalsnpVcf[@snpVcf] == 2)
			{
				 $v1  = 1;
			}
			if(rand() < 0.5 && $v1 ==0)
			{
				#print "test1 $temp2\n";
				@temp = split('',$temp2);
				$beforeinv = "";
				$inv="";
				$afterinv = "";
				for($i=0;$i<@temp;$i++)
				{
					if($i<=$eventstart)
					{
							$beforeinv = $beforeinv.$temp[$i];
					}
					elsif($i>$eventstart && $i<=$eventstop)
					{
						$inv = $inv.$temp[$i];
					}
					else
					{
						$afterinv = $afterinv.$temp[$i];
					}
					
				}
				$inv = reverse $inv;
				$inv =~ tr/ACGTacgt/TGCAtgca/;
				$temp2 = $beforeinv.$inv.$afterinv;
				@temp = split('',$temp2);
				for($i=0;$i<@temp;$i++)
				{
					if((($i+1)%$fasta_len) == 0)
						{
							#print "success\n";
							$temp[$i] = $temp[$i]."\n";
						}
				}
				$temp2 = join('',@temp);
				$temp2 =~ s/\n$//g;
				#print "test2 $temp2\n";
				$line=$temp2;
				push(@snpVcf,1);
				#print @snpVcf."\n";
			}	
			else
			{
				 #print @snpVcf."\n";
				#print "$v1\t$flag"."\t".$maternalsnpVcf[@snpVcf]."\n";
				if($v1 !=0)
				{
					print LOG "ignore event $$_[0] multiple events occuring on the same lines of $fasta_len bp length\n";
					$temp_fail =~ s/\n$//g;
					$line=$temp_fail;	
					push(@snpVcf,2);
				}
				elsif($flag ==1 && $maternalsnpVcf[@snpVcf] == 0)
				{
					#print "sucess $$_[0] $maternalsnpVcf[@snpVcf]\n";
					@temp = split('',$temp2);
					$beforeinv = "";
					$inv="";
					$afterinv = "";
					for($i=0;$i<@temp;$i++)
					{
						if($i<=$eventstart)
						{
								$beforeinv = $beforeinv.$temp[$i];
						}
						elsif($i>$eventstart && $i<=$eventstop)
						{
							$inv = $inv.$temp[$i];
						}
						else
						{
							$afterinv = $afterinv.$temp[$i];
						}
					
					}
					$inv = reverse $inv;
					$inv =~ tr/ACGTacgt/TGCAtgca/;
					$temp2 = $beforeinv.$inv.$afterinv;
					@temp = split('',$temp2);
					for($i=0;$i<@temp;$i++)
					{
						if((($i+1)%$fasta_len) == 0)
						{
							#print "success\n";
							$temp[$i] = $temp[$i]."\n";
						}
					}
					$temp2 = join('',@temp);
					$temp2 =~ s/\n$//g;
					$line=$temp2;
					push(@snpVcf,1);
				}  
				else
				{	
					$temp_fail =~ s/\n$//g;
					$line=$temp_fail;
					push(@snpVcf,0);
				}
			}

	
	$temp = $$_[1];		
	$index++;
}

if($temp ne "")
{
	print $maternal $line."\n";
        #print $paternal $line."\n";
        while($line=<$buff>)
        {
        	print $maternal $line;
            #print $paternal $line;
        }
}

close($buff);
close($maternal);
return(@snpVcf);

}

############################################TRANSLOCATIONS##########################################
sub trans{


@maintrans=();
@trans1 = ();
@trans2 = ();
$transFlag_doulerc = 0;
#@trans3 = ();
#@trans4 = ();

$temp = "";
$transindex = -1;
my($val) = shift @_;
@snpVcf=();
#looping through translocation array
	for($index=0;$index<@translocationsHash;$index++)
	{
	
		#trans location separate index as two events considered as one
		$transindex++;
		$v1 = 0;
		#transflag to indicate the reciproal translocation as 4 events represent as one
		if($transFlag != 1 && $transFlag != 4)
		{	
			$transFlag =0;
		}
		#trans 1 and trans2 array represents two events as one
		@trans1 = @{$translocationsHash[$index]};
		$index++;
		@trans2 = @{$translocationsHash[$index]};
		#print "$trans1[4]\n$trans2[4]\n";
		
		#extract second chromosome and position (target)
		$Transaltchr = $trans1[4];
		$Transaltchr =~ s/\:.*//g;
		@Transaltchr = ($Transaltchr =~ m/(\d+)/g);
		#print $trans1[4]."\t$trans2[4]\n";
		#seleting the main event based on the type (selecting one from two events to implement)
		if(($trans1[4] =~ m/^[ATGCatgc]/ && $trans2[4] =~ m/^[ATGCatgc]/) || ($trans1[4] =~ m/^[ATGCatgc]/ && $trans2[4] !~ m/^[ATGCatgc]/))
		{
			@maintrans = @trans1;
		}
		elsif($trans1[4] !~ m/^[ATGCatgc]/ && $trans2[4] =~ m/^[ATGCatgc]/)
		{
			@maintrans = @trans2;
		}
		elsif($trans1[4] =~ m/[ATGCatgc]$/ && $trans2[4] =~ m/[ATGCatgc]$/)
		{
			@maintrans = @trans1;
			$transFlag_doulerc = 1;
			#die;
		}
		else
		{
			print LOG "ignoring events $trans1[0] && $trans2[0] unexpected translocation\n";
		}
		
		#if first reciprocal event failed, second one shud fail
		if($transFlag == 4)
		{
			@maintrans =();
			print LOG "ignoring event $trans1[0] reciprocal translocation as first one ignored\n";
			$transFlag = 5;
		}
		
		#checking if the event occured on the same chr before
		if($transFlag == 0 && @maintrans > 0 && (exists($translocationsTrackHash{$trans1[1]}) || exists($translocationsTrackHash{$Transaltchr[0]})))
		{
			print LOG "ignoring event $trans1[0] multiple translocations cannot occur on same chromosome\n";
			@maintrans =();
			#print $trans1[5]."\n";
			if($trans1[5]  =~ m/reciprocal/ &&  $transFlag == 0)
			{
				$transFlag = 4;
			}
		}
		else
		{
			$translocationsTrackHash{$trans1[1]} =1;
			$translocationsTrackHash{$Transaltchr[0]} =1;
			
		}
		
		#reciprocal event handle
		if($trans1[5]  =~ m/reciprocal/ &&  $transFlag != 4)
		{
			$transFlag++;
		}
		#if(exists($translocationsTrackHash{$maintrans
		#print "final test trans $index $maternalsnpVcf[$index]\n";
		
		#event enter main loop after pass all the filter steps
		if(@maintrans > 0 )
		{
			#print "TransFlag:$transFlag\n";
		

			#print "@maintrans\n";
			$k = 0;
			#print "final test trans $index $maternalsnpVcf[$index]\n";
			
			#if maternal event fails an event becoz of random probability atleast it shud occur in parental chr
			if($flag ==1 && $maternalsnpVcf[$transindex]==0)
			{
				$k =1;
			}
			$rand = rand();
			#print "$transFlag\t$transindex\n";
			
			#reciprocal event handling if first passes second reciprocal shud also pass or if first event fails second shud fail
			if($transFlag == 2 && $snpVcf[$transindex-1] ==1)
			{
				$k =1;
			}
			elsif($transFlag == 2 && $snpVcf[$transindex-1] ==0)
			{
				$rand = 0.51;
			}
			else
			{
			
			}
			# if event passes the random probability condition implements the event
			if($rand <0.5 || $k == 1)
			{
				
					$buff = "CHR$maintrans[1]";
					$maternal = "MATERNAL$maintrans[1]";
				if($transFlag_doulerc == 0)
				{
					#selecting the previous event file to read
					if(exists($snpsOutFilenames{$maintrans[1]}))
					{
					
						$file = "$tempdir/$val"."_$snpsOutFilenames{$maintrans[1]}";
					}
					else
					{
						$file = "$refgenome_dir/chr$maintrans[1]".'.fa';
					}
					#print "$file\n";
					open $buff,"$file" or die "$file\n";
				
					#$file = "TRANSchr$maintrans[1]".'.fa';
					#if($flag== 1)
					#{
					#	$snpsOutFilenames{$maintrans[1]} =  $file;
					#}
				#if event is double reciprocal or other 
				
					#creating the new translocation file to write
					$file = "$val"."_TRANSchr$maintrans[1]".'.fa';
					open($maternal,">$tempdir/$file");
					$line=<$buff>;
					print $maternal $line;
					$length =0;
					$prevlength =0;
					
					#creating left over chr for non reciprocal translocations
					if($transFlag == 0)
					{
							open(LEFTOVER,">$tempdir/leftover_$file");
							print LEFTOVER ">chr$maintrans[1]\n";
					}
					
					#looping through the file to reach the position
					while($maintrans[2] > $length)
					{
				
						if($length != 0)
						{
							print $maternal $line."\n";
						}
						if($line=<$buff>)
						{
			
							chomp($line);
							$prevlength = $length;
							#$length=$length+length($line);
							$length=$length+$fasta_len;
						}
						else
						{
							die "event $maintrans[0] coord exceeds chr length please change input vcf file and rerun ignoring event not possible (End of chr file)\n";
						}
					}	
					#print "test test $file $line\n";
					
					#seleting the exact position
					$temp  = $maintrans[2]%$fasta_len;
					#print $temp."\n";
					@temp = split('',$line);
					
					
					
					for($i=0;$i<@temp;$i++)
					{
						if($i<$temp)
						{
							print $maternal $temp[$i];
						}
						elsif($transFlag == 0)
						{
							print LEFTOVER $temp[$i];
						}
						else
						{
						
						}
					}
					print $maternal "\n";
					
					#writing remaining lines to leftover chr
					if($transFlag == 0)
					{
						print LEFTOVER "\n";
						while($line=<$buff>)
						{
							print  LEFTOVER $line;
						}
					}
					close($buff);
					
					#based on the type of translocation decided need to flip the sequence or not
					if($maintrans[4] =~ m/\[\d+:\d+\[/)
					{
						#print "category1\n";
						
						#selecting the target chr and position
						$Transaltchr = $maintrans[4];
						#$Transaltchr =~ s/\:.*//g;
						@Transaltchr = ($Transaltchr =~ m/(\d+)/g);
						
						#selecting the target chr file
						if(exists($snpsOutFilenames{$Transaltchr[0]}))
						{
					
							$file = "$tempdir/$val"."_$snpsOutFilenames{$Transaltchr[0]}";
						}
						else
						{
							$file = "$refgenome_dir/chr$Transaltchr[0]".'.fa';
						}
						#print LOG "maintestTrans\t".$file."\t$Transaltchr[1]\n";
						
						#reading fron the target file
						open BUFF,"$file" or die "no file exists $file\n";
						$line = <BUFF>;
						
						$file = "$val"."_TRANSchr$Transaltchr[0]".'.fa';
						
						#creating the left over target chr if the event is non reciprocal event
						if($transFlag == 0)
						{	
							open(LEFTOVER,">$tempdir/leftover_$file");
							print LEFTOVER ">chr$Transaltchr[0]\n";
						}
						
						#selecting the exact target position
						$temp1 = int(($Transaltchr[1]-0.1)/$fasta_len);
						$temp  = $Transaltchr[1]%$fasta_len;
						#print "final final $temp1 $temp\n";
						$num =1;
						while ($line = <BUFF>)
						{
							#print "omit $line" ;
							#for non reciprocal event writing the left over chr
							if($transFlag == 0)
							{	
								print LEFTOVER $line;
							}
							if($num == $temp1+1)
							{
								last;
							}
							
							$num++;
						}
						chomp($line);
						@temp = split("",$line);
						for($i=0;$i<@temp;$i++)
						{
							if($i>$temp-1)
							{
								print $maternal $temp[$i];
							}
							#for non reciprocal event writing the left over chr
							elsif($transFlag == 0)
							{
								print LEFTOVER $temp[$i];
							}
							else
							{
							
							}
						}
						print $maternal "\n";
						if($transFlag == 0)
						{
							print LEFTOVER "\n";
							close(LEFTOVER);
						}
						while ($line = <BUFF>)
						{
							print $maternal "$line" ;
							
						}
						
						close(BUFF);
						close($maternal);
						
					}
					
					#if the event is second type of reciprocal event then the target shud be fliped
					else
					{	
					
						#selecting the target chr and position
						$Transaltchr = $maintrans[4];
						#$Transaltchr =~ s/\:.*//g;
						@Transaltchr = ($Transaltchr =~ m/(\d+)/g);
						#print $file."\n";
						#print "category2\n";
						
						#sleecting the target file to read
						if(exists($snpsOutFilenames{$Transaltchr[0]}))
						{
					
							$file = "$tempdir/$val"."_$snpsOutFilenames{$Transaltchr[0]}";
						}
						else
						{
							$file = "$refgenome_dir/chr$Transaltchr[0]".'.fa';
						}
						#print $file."\t$Transaltchr[1]\n";
						#open TRANS 
						
						#reading the target file backwords
						$fh = File::ReadBackwards->new($file)  or die "can't read $file: $!\n";
						$num =-1;
						#counting the number of lines in the file
						while ( defined($line = $fh->readline) )
						{
							#print "my $line" ;
							$num++;
						}
						
						#selecting the exact line and position
						$temp1 = $num-int(($Transaltchr[1]-0.1)/$fasta_len);
						$temp  = $fasta_len-$Transaltchr[1]%$fasta_len;
						#print $temp."\t$temp1\t$num\n";
						
						#readin the file again to select the exact line
						$fh = File::ReadBackwards->new($file)  or die "can't read $file: $!\n";
						$num1 = 0;
						while ( defined($line = $fh->readline) )
						{
							if($num1 == $temp1-1)
							{
								last;
							}
							$num1++;
						}		
						#print "final $line\n";
						chomp($line);
						
						#select exact line and pisition the do reverse complement and print in the output file
						@temp = split('',$line);
						$str= "";
						for($i=$temp;$i<@temp;$i++)
						{
							$str=$str.$temp[$i];
						}
						$str = reverse($str);
						$str =~ tr/ACGTacgt/TGCAtgca/;
						print $maternal "$str\n";
						
						$num1 = $num-$num1;
						#print "num $num num1 $num1\n";
						while ( defined($line = $fh->readline) )
						{
							
							if($num1 == 1)
							{
								last;
							}
							chomp($line);
							$str = $line;
							$str = reverse($str);
							$str =~ tr/ACGTacgt/TGCAtgca/;
							print $maternal "$str\n";
							$num1--;
							#close($maternal);
							
						}
						close($maternal);
						
						#left over part of the second chromosome
						if($transFlag == 0)
						{
						
							open BUFF,"$file" or die "no file exists $file\n";
							$line = <BUFF>;
							$file = "$val"."_TRANSchr$Transaltchr[0]".'.fa';
							open(LEFTOVER,">$tempdir/leftover_$file");
							print LEFTOVER ">chr$Transaltchr[0]\n";
							$temp1 = int(($Transaltchr[1]-0.1)/$fasta_len);
							$temp  = $Transaltchr[1]%$fasta_len;
							#print "final final $temp1 $temp\n";
							$num =1;
							while ($line = <BUFF>)
							{
								#print "omit $line" ;
								if($num == $temp1)
								{
									last;
								}
								$num++;
							}
							chomp($line);
							@temp = split("",$line);
							for($i=0;$i<@temp;$i++)
							{
								if($i>$temp-1)
								{
									print LEFTOVER $temp[$i];
								}
							}
							print LEFTOVER "\n";
							while ($line = <BUFF>)
							{
								print LEFTOVER "$line" ;
								
							}
							close(LEFTOVER);
							close(BUFF);
						}	
					}
					push(@snpVcf,1);
				}
				else
				{
						$Transaltchr = $maintrans[4];
						@Transaltchr = ($Transaltchr =~ m/(\d+)/g);	
						$file = "$val"."_TRANSchr$maintrans[1]".'.fa';
						#print "writefile $file\n";
						open($maternal,">$tempdir/$file");
						print $maternal ">chr$maintrans[1]\n";
						if(exists($snpsOutFilenames{$Transaltchr[0]}))
						{
					
							$file = "$tempdir/$val"."_$snpsOutFilenames{$Transaltchr[0]}";
						}
						else
						{
							$file = "$refgenome_dir/chr$Transaltchr[0]".'.fa';
						}
						
						$fh = File::ReadBackwards->new($file)  or die "can't read $file: $!\n";
						$num =-1;
						#counting the number of lines in the file
						while ( defined($line = $fh->readline) )
						{
							$num++;
						}
						
						#selecting the exact line and position
						$temp1 = $num-int(($Transaltchr[1]-0.1)/$fasta_len);
						$temp  = $fasta_len-$Transaltchr[1]%$fasta_len;
						#print "final test main $temp1 $temp\n";
						#die;
						#readin the file again to select the exact line
						$fh = File::ReadBackwards->new($file)  or die "can't read $file: $!\n";
						$num1 = 0;
						while ( defined($line = $fh->readline) )
						{
							if($num1 == $temp1-1)
							{
								last;
							}
							chomp($line);
							$line = reverse($line);
							$line =~ tr/ACGTacgt/TGCAtgca/;
							print $maternal $line."\n";
							$num1++;
						}	
						chomp($line);
						@temp = split("",$line);
						$temp5 = "";
						for($i=0;$i<@temp;$i++)
						{
							if($i<$temp)
							{
								$temp5 = $temp5.$temp[$i];
							}
						}
						#print "hi $temp5\t$line\t$temp\n";
						$temp5 = reverse($temp5);
						$temp5 =~ tr/ACGTacgt/TGCAtgca/;
						print $maternal $temp5."\n";
						#print $maternal "end\n";
						#creating leftover chr
						$file = "$val"."_TRANSchr$Transaltchr[0]".'.fa';
						open(LEFTOVER,">$tempdir/leftover_$file");
						#print "LEFTOVER $tempdir/leftover_$file\n";
						#print LEFTOVER ">chr$Transaltchr[0]\n";
						$temp1 = int(($Transaltchr[1]-0.1)/$fasta_len);
						$temp  = $Transaltchr[1]%$fasta_len;
						#print "final test main $temp1 $temp\n";
						if(exists($snpsOutFilenames{$Transaltchr[0]}))
						{
					
							$file = "$tempdir/$val"."_$snpsOutFilenames{$Transaltchr[0]}";
						}
						else
						{
							$file = "$refgenome_dir/chr$Transaltchr[0]".'.fa';
						}
						#print "File read $file\n";
						open(BUFF,$file) or die "no file exists $file\n";
						$line = <BUFF>;
						print LEFTOVER $line;
						$num1=0;
						while ($line=<BUFF>)
						{
							if($num1 == $temp1)
							{
								last;
							}
							chomp($line);
							print LEFTOVER $line."\n";
							$num1++;
						}	
						chomp($line);
						@temp = split("",$line);
						$temp5 = "";
						for($i=0;$i<@temp;$i++)
						{
							if($i<$temp)
							{
								$temp5 = $temp5.$temp[$i];
							}
						}
						print LEFTOVER $temp5."\n";
						close(LEFTOVER);
						close(BUFF);
						if(exists($snpsOutFilenames{$maintrans[1]}))
						{
					
							$file = "$tempdir/$val"."_$snpsOutFilenames{$maintrans[1]}";
						}
						else
						{
							$file = "$refgenome_dir/chr$maintrans[1]".'.fa';
						}
						open(BUFF,$file) or die "no file exists $file\n";
						#creating leftover chr
						$file = "$val"."_TRANSchr$maintrans[1]".'.fa';
						open(LEFTOVER,">$tempdir/leftover_$file");
						$line=<BUFF>;
						print LEFTOVER $line;
						$temp1 = int(($maintrans[2]-0.1)/$fasta_len);
						$temp  = $maintrans[2]%$fasta_len;
						$num1=0;
						while ($line=<BUFF>)
						{
							if($num1 == $temp1)
							{
								last;
							}
							chomp($line);
							print LEFTOVER $line."\n";
							$num1++;
						}
						chomp($line);
						@temp = split("",$line);
						for($i=0;$i<@temp;$i++)
						{
							if($i<$temp)
							{
								#$temp5 = $temp5.$temp[$i];
								print LEFTOVER $temp[$i];
							}
							else
							{
								print $maternal $temp[$i];
							}
						}
						print LEFTOVER "\n";
						print $maternal "\n";
						while($line=<BUFF>)
						{
							print $maternal $line;
						}
						close($maternal);
						close($buff);
						close($maternal);
						close(LEFTOVER);
						push(@snpVcf,1);
				}
			}
			#if randome probability is less than 0.50 the push 1 in to array
			else
			{
				push(@snpVcf,0);
			}
		}
		else
		{
			push(@snpVcf,2);
		}
	}
	return(@snpVcf);
}

close(LOG);
