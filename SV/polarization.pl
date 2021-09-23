#!/usr/bin/perl
use strict;
use warnings;

my $p2pdels = $ARGV[0]; #  pseIns_per2pse.bed
my $p2pins = $ARGV[1]; # perDel_per2pse.bed


print "block\tDpseDelChrom\tDelstart\tDelend\tsvimScore\tDperInsChrom\tInsstart\tInsend\tsvimScore\tDelinsvmu Insinsvmu\tDelinsvmu>Dmir Insinsvmu>Dmir\tDelinsvmu=!Dmir Insinsvmu=!Dmir\tDelnotinsvmu Insnotinsvmu\n";
#print "block\tDpseInsChrom\tInsstart\tInsend\tsvimScore\tDperDelChrom\tDelstart\tDelend\tsvimScore\tInsinsvmu Delinsvmu\tInsinsvmu>Dmir Delinsvmu>Dmir\tInsinsvmu=!Dmir Delinsvmu=!Dmir\tInsnotinsvmu Delnotinsvmu\n";

my %var; my $indmir;
open FILE1, "$p2pdels" || die "Can't open FILE1\n";
while(my $line = <FILE1>){
	chomp($line);
	my @campos = split("\t",$line);
	my @campos2 = split(" ",$campos[3]);
	my $block = $campos2[2];
	my $blockrange = $campos2[1] * .25;
	my $minr = $campos2[1] - $blockrange;
	my $maxr = $campos2[1] + $blockrange;
	my $psedel = $campos[0] . "\t" . $campos[1] . "\t" . $campos[2] . "\t" . $campos2[0];
	my $psecol = $campos[0].$campos[1].$campos[2];
	#print "$psecol\t";
	my $perin;
	my $percol;
	my @perins = `awk '{ if (\$NF == "$block") print \$0 }' $p2pins`; chomp(@perins);
	if (scalar @perins > 1) {
		foreach my $val(@perins){
			my @percam = split(" ", $val);
			if ($percam[4] > $minr && $percam[4] < $maxr){
				#print "$val\n";
				$perin = $percam[0] . "\t" . $percam[1] . "\t" . $percam[2] . "\t" . $percam[3];
				$percol = $percam[0].$percam[1].$percam[2];
				#print "$percol\n";
			}
		}
	}
	if (scalar @perins == 1){
		my @percam = split(" ", $perins[0]);
		#print "$perins[0]\n";
		$perin = $percam[0] . "\t" . $percam[1] . "\t" . $percam[2] . "\t" . $percam[3];
		$percol = $percam[0].$percam[1].$percam[2];
		#print "$percol\n";
	}

	#my $fullblock = $block . "\t" . $psedel ."\t". $perin;


	############## Polarize Dpse Dels  in Dmir ####################################################################################
	#print "$fullblock\n";
	my $dpsedelindmir = "no";
	my @dpseconfdmir = `awk '{ var = \$1\$2\$3; if (var == "$psecol") print \$0 }' Transformed_pseInsper2pse_svmuconfirmed_R.txt`;chomp(@dpseconfdmir);
	my $dpsedelindmirdiffsize = "no";
	my @dpseconfdmirdf = `awk '{ var = \$1\$2\$3; if (var == "$psecol") print \$0 }' Transformed_pseInsper2pse_svmuconfirmed_inDmir.txt`;chomp(@dpseconfdmirdf);
	my $dpsecollinearblockindmir = "no";
	my @dpsecoldf = `awk '{ var = \$1\$2\$3; if (var == "$psecol") print \$0 }' Transformed_pseIns_svmuper2pse.txt_notsvmu_confirmed.txt`;chomp(@dpsecoldf);
	my $dpsenotfoundindmir = "no";
	my @dpsenotfoundmir = `awk '{ var = \$1\$2\$3; if (var == "$psecol") print \$0 }'  Transformed_pseIns_per2pse.bed_notsvmu_confirmed.txt`;chomp(@dpsenotfoundmir);

	if (scalar @dpseconfdmir == 1){
		$dpsedelindmir = "yes";
	}
	if (scalar @dpseconfdmirdf == 1){
		$dpsedelindmirdiffsize = "yes";
	}
	if (scalar @dpsecoldf == 1){
		$dpsecollinearblockindmir = "yes";
	}
	if (scalar @dpsenotfoundmir == 1){
		$dpsenotfoundindmir = "yes";
	}


	my $dperinsindmir = "no";
	my @dperconfdmir = `awk '{ var = \$1\$2\$3; if (var == "$percol") print \$0 }' Transformed_perDelper2pse_svmuconfirmed_R.txt`;chomp(@dperconfdmir);
	my $dperinsindmirdiffsize = "no";
	my @dperconfdmirdf = `awk '{ var = \$1\$2\$3; if (var == "$percol") print \$0 }' Transformed_perDelper2pse_svmuconfirmed_inDmir.txt`;chomp(@dperconfdmirdf);
	my $dpercollinearblockindmir = "no";
	my @dpercoldf = `awk '{ var = \$1\$2\$3; if (var == "$percol") print \$0 }' Transformed_perDel_svmuper2pse.txt_notsvmu_confirmed.txt`;chomp(@dpercoldf);
	my $dpernotfoundindmir = "no";
	my @dpernotfoundmir = `awk '{ var = \$1\$2\$3; if (var == "$percol") print \$0 }' Transformed_perDel_per2pse.bed_notsvmu_confirmed.txt`;chomp(@dpernotfoundmir);

	if (scalar @dperconfdmir == 1){
		$dperinsindmir = "yes";
	}
	if (scalar @dperconfdmirdf == 1){
		$dperinsindmirdiffsize = "yes";
	}
	if (scalar @dpercoldf == 1){
		$dpercollinearblockindmir = "yes";
	}
	if (scalar @dpernotfoundmir == 1){
		$dpernotfoundindmir = "yes";
	}

	my $fullblock = $block ."\t". $psedel ."\t". $perin ."\t". $dpsedelindmir ." ". $dperinsindmir ."\t". $dpsedelindmirdiffsize ." ". $dperinsindmirdiffsize ."\t". $dpsecollinearblockindmir ." ". $dpercollinearblockindmir ."\t". $dpsenotfoundindmir ." ". $dpernotfoundindmir;
	print "$fullblock\n";

}
close FILE1;

#print "$var{block_1}\n";
