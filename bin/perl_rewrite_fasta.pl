$file = $ARGV[0];
chomp($file);
$fasta_len = 60;
open(BUFF,$file) or die " no file found\n";
$line=<BUFF>;
print $line;
$temp1 = 0;
while($temp=<BUFF>)
{
	chomp($temp);
	$temp =~ s/ //g;
	#if($temp ne "")
	#{
	#	print WRBUFF "$temp\n";
	#}
	
	undef(@array);
	@array = split('',$temp);
	for($i=0;$i<@array;$i++)
	{
		if($temp1 == $fasta_len)
		{
			print  "\n";
			$temp1 = 0;
		}
		print $array[$i];
		$temp1++;
	}
}
if($temp1 != 0)
{
	print "\n";
}
