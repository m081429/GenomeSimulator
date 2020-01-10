$nucleo[$inLinePos]="C";
print randInsertionArray($nucleo[$inLinePos])."\n";
sub randInsertionArray{
	my $max_base = 50;
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
sub random_int_between {
		my($min, $max) = @_;
		# Assumes that the two arguments are integers themselves!
		return $min if $min == $max;
		($min, $max) = ($max, $min)  if  $min > $max;
		return $min + int rand(1 + $max - $min);
		}
=head
print random_int_between(10,20)."\n";
sub random_int_between 
{
	my($min, $max) = @_;
	# Assumes that the two arguments are integers themselves!
	return $min if $min == $max;
	($min, $max) = ($max, $min)  if  $min > $max;
	return $min + int rand(1 + $max - $min);
}
