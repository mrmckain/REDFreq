#!/usr/bin/perl 
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);

my $help = 0;
my $bin_size = 50;
my $genome_file;
my $resenz;

GetOptions('help|?' => \$help, "genome=s" => \$genome_file, "enzymes=s" => \$resenz, "bin=i" => \$bin_size);

print "Options:\n --genome = Sequence file to estimate restiction digest fragment sizes.\n--enzymes = Comma delimited list of enzymes to use.  Enzymes should be in the format: EcoRI, where capitalization matters. If you are unsure, check the file TypeII_Restiction_Enzymes.txt.\n--bin = Size of bins. Default is 50bp.\n" and exit if $help;


my %alt_nuc;
push(@{$alt_nuc{K}}, "G");
push(@{$alt_nuc{K}}, "T");
push(@{$alt_nuc{M}}, "A");
push(@{$alt_nuc{M}}, "C");
push(@{$alt_nuc{R}}, "A");
push(@{$alt_nuc{R}}, "G");
push(@{$alt_nuc{Y}}, "C");
push(@{$alt_nuc{Y}}, "T");
push(@{$alt_nuc{S}}, "C");
push(@{$alt_nuc{S}}, "G");
push(@{$alt_nuc{W}}, "A");
push(@{$alt_nuc{W}}, "T");
push(@{$alt_nuc{B}}, "C");
push(@{$alt_nuc{B}}, "G");
push(@{$alt_nuc{B}}, "T");
push(@{$alt_nuc{V}}, "C");
push(@{$alt_nuc{V}}, "G");
push(@{$alt_nuc{V}}, "A");
push(@{$alt_nuc{H}}, "C");
push(@{$alt_nuc{H}}, "A");
push(@{$alt_nuc{H}}, "T");
push(@{$alt_nuc{D}}, "A");
push(@{$alt_nuc{D}}, "G");
push(@{$alt_nuc{D}}, "T");
push(@{$alt_nuc{D}}, "A");
push(@{$alt_nuc{D}}, "T");
push(@{$alt_nuc{D}}, "G");
push(@{$alt_nuc{D}}, "C");

my @enz = split(/,/, $resenz);
my %enzymes = map {$_ => 1} @enz;
my %sites;

open my $enz_files, "<", $Bin . "/TypeII_Restiction_Enzymes.txt";
while(<$enz_files>){
		chomp;
		my @tarray = split /\s+/;
		if(exists $enzymes{$tarray[1]}){
			my @possible_sites = site_generator($tarray[0]);
			for my $psite (@possible_sites){
				push(@{$sites{$tarray[1]}}, $psite);
				my $tempsite = $psite;
				$tempsite =~ tr/ATCGatcg/TAGCtagc/;
				$tempsite = reverse($tempsite);
				push(@{$sites{$tarray[1]}}, $tempsite);
			}
		}
}



my %genome;
my $id;
open my $IN, "<", $genome_file or die "Cannot open $genome_file\n";
while(<$IN>){
	chomp;
	if(/>/){
		$id=$_;
	}
	else{
		$genome{$id}.=$_;
	}
}



my %positions;
for my $ids (sort keys %genome){
	my $offset = 0;
	for my $enzyme (sort keys %sites){
		for my $cut (@{$sites{$enzyme}}){
			$offset = 0;
			#print "cut is $cut\n";
			while($offset<length($genome{$ids})){
				my $pos = index($genome{$ids}, $cut, $offset);
			#	print "$pos\n";
				if($pos==-1){
					$offset= length($genome{$ids});
				}
				unless($pos==-1){
					$positions{$ids}{$pos}=$enzyme;
					$offset = $pos +1;
				}
			}
		}
	}
}

%genome=();
my @cutsizes;
my $prev;
my $prev_enzyme;
for my $idz (sort keys %positions){
	$prev=-1;
	for my $val (sort {$a<=>$b} keys %{$positions{$idz}}){
		#print "val is $val\n";
		if(scalar keys %enzymes == 1){
			 my $diff = $val-$prev;
                #       print "$diff\n";
                        push(@cutsizes, $diff);
                        $prev=$val;
                        $prev_enzyme=$positions{$idz}{$val};
                }
		else{
		if($prev >= 0 && $prev_enzyme ne $positions{$idz}{$val}){
			my $diff = $val-$prev;
		#	print "$diff\n";
			push(@cutsizes, $diff);
			$prev=$val;
			$prev_enzyme=$positions{$idz}{$val};
		}
		else{
			$prev=$val;
			$prev_enzyme=$positions{$idz}{$val};
		}}
	}
}

@cutsizes = sort {$a <=> $b} @cutsizes;
my $max = $cutsizes[$#cutsizes];
my $lmax = int($max/$bin_size);
my $bmax = ($lmax+1) *$bin_size;
my %frequencies;
for(my $i=$bin_size; $i <= $bmax; $i+=$bin_size){
	$frequencies{$i}=0;
}
for(my $i=0; $i <= $bmax; $i+=$bin_size){
	for my $cusize (@cutsizes){
		if($cusize <= $i && $cusize >= $i-$bin_size){
			$frequencies{$i}++;
		}
	}
}

my $enz = join(".", @enz);
open my $OUT, ">", "$genome_file" . ".$enz" . ".$bin_size\_freqs.txt";
for my $freqs (sort {$a<=>$b} keys %frequencies){
	print $OUT "$freqs\t$frequencies{$freqs}\n";
}
close $OUT;

open my $OUT2, ">", "$genome_file" . ".$enz" . "_Summary.txt";
for my $ids (sort keys %positions){
	my $count = scalar keys %{$positions{$ids}};
	print $OUT2 "$ids\t$count\n";
}

###########################
sub site_generator {
	
	my $seq = $_[0];

	my %temp_seqs;

	my @split_seq = split(//, $seq);
	for my $base (@split_seq){
			if(exists $alt_nuc{$base}){
					if(%temp_seqs){
							for my $tseq (keys %temp_seqs){
								for my $nbase (@{$alt_nuc{$base}}){
									my $nseq .= $tseq . $nbase;
									$temp_seqs{$nseq}=1;
								}
								delete $temp_seqs{$tseq};
							}
					}
					else{
						for my $nbase (@{$alt_nuc{$base}}){
							$temp_seqs{$nbase}=1;
						}
					}
			}
			else{
				if(%temp_seqs){
						for my $tseq (keys %temp_seqs){
							my $nseq .= $tseq . $base;
							delete $temp_seqs{$tseq};
							$temp_seqs{$nseq}=1;
						}
				}
				else{
					$temp_seqs{$base}=1;
				}

			}
	}
	my @final_seqs;
	for my $seqs (keys %temp_seqs){
		push(@final_seqs, $seqs);
	}
	return @final_seqs;
}


