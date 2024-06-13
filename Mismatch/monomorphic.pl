
#!/usr/bin/perl

use strict;
use warnings;

my $Tab = shift;

if ($Tab =~ /.gz$/) {

	open(IN, "gunzip -c $Tab |") || die "can’t open pipe to $Tab";

}
else {

	open(IN, $Tab) || die "can’t open $Tab";

}

<IN>;

while(<IN>){

	if($. % 1000000 == 0){print STDERR $.."\n";}

	my @tt = split /\t/;

	shift @tt;

	my $key = shift @tt;

	shift @tt;

	shift @tt;

	@tt = grep !/0 0/, @tt;

	if (keys %{{ map {$_, 1} @tt }} == 1) {

		#print STDERR $.."\n";

		#print $key."\n";

		my @tt2 = split / /, $tt[0];

		if($tt2[0] eq $tt2[1]){

			print STDERR $.."\t".$key."\n";

			print $key."\n";

		}

	}

}

close IN;

