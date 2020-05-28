#!/bin/usr/perl -w
use warnings;
use Getopt::Long;

my %opt = ('ksize' => 4);

GetOptions(\%opt, 'AsOne|a', 'ksize|k:i', 'canonical_kmer|c');

@ARGV == 1 || die " Script for kmer counting.
	perl $0 seq.fna [options] \n
	--AsOne            -a  if multi-sequences in a fasta file, count as one	\n
	--ksize            -k  k-mer size, default is 4\n
	--canonical_kmer   -c  canonical k-mer type\n";

my $fasta = shift;

open FA, $fasta || die $!;

$/ = '>'; <FA>;

my(%kmer_count, %seqID2TNF, @seqIDs);

while(<FA>){
	chomp;
	my @temp = split /\n/;
	my $id = shift @temp;
	my $seq = uc(join("", @temp));
	my $len = length($seq);
	if($opt{AsOne}){
		for my $i (0..$len-$opt{ksize}){
			my $kmer = substr($seq, $i, $opt{ksize});
			$kmer =~ /N/ && next;
			$kmer_count{$kmer}++;
		}
	}else{
		my %temp_tnf = kmer_frequence($seq, $opt{ksize});
		my %cano_tnf = $opt{canonical_kmer} ? com_rev_kmer_num_add(\%temp_tnf) : %temp_tnf;
		for my $k (keys %cano_tnf){
			$seqID2TNF{$k}->{$id} = $cano_tnf{$k};
		}
		push @seqIDs, $id;
	}
}
$/ = "\n";
close FA;

if($opt{AsOne}){
	if(!$opt{canonical_kmer}){
		for my $k (sort { $a cmp $b } keys %kmer_count){
			print "$k\t$kmer_count{$k}\n";
		}
	}else{
		my %nkmer_count = com_rev_kmer_num_add(\%kmer_count);
		for my $k (sort { $a cmp $b } keys %nkmer_count){
			print "$k\t$nkmer_count{$k}\n";
		}
	}
}else{
	my @kmers = sort { $a cmp $b } keys %seqID2TNF;
	print join("\t", "#ID", @seqIDs), "\n";
	for my $k (@kmers){
		my %inner_hash = %{$seqID2TNF{$k}};
		my @contigs = sort { $a cmp $b } keys %inner_hash;
		my @temp_vals;
		for my $id (@seqIDs){
			my $temp_val = $inner_hash{$id} ? $inner_hash{$id}  : 0;
			push @temp_vals, $temp_val;
		}
		print $k, "\t", join("\t", @temp_vals), "\n";
	}
}



sub kmer_frequence{
	my ($seq, $ksize) = @_;
	$seq = uc($seq);
	my $len = length($seq);
	my %count_km;
	for my $i (0..$len-$ksize){
		my $kmer = substr($seq, $i, $ksize);
		$kmer =~ /N/ && next;
		$count_km{$kmer}++;
	}
	return %count_km;
}

sub com_rev{
	my ($seq) = @_;
	my $com_seq = $seq;
	$com_seq =~tr/ATGC/TACG/;
	my $cr_seq = reverse($com_seq);
	return $cr_seq;
}


sub com_rev_kmer_num_add{
	my ($kmer_num) = @_;
	my %alpha_order = ('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4);
	my %kmer_exist; my %nkmer_count;
	my %kmer_count_hash = %{$kmer_num};
	for my $k (sort { $a cmp $b } keys %kmer_count_hash){
		my $crk = com_rev($k);
		$n = ($kmer_count_hash{$crk} && $k ne $crk) ? $kmer_count_hash{$crk} : 0;
		if(!$kmer_exist{$k}){
			$kmer_exist{$k} = 1;
			$kmer_exist{$crk} = 1;
		}else{
			next;
		}
		my $k_first_letter = substr($k, 0, 1);
		my $crk_first_letter = substr($crk, 0, 1);
		#print $k_first_letter, "\t", $crk_first_letter, "\n";
		my $adjk = $alpha_order{$k_first_letter} <= $alpha_order{$crk_first_letter} ? $k : $crk;
		$nkmer_count{$adjk} = $kmer_count_hash{$k} + $n;
	}
	return %nkmer_count;
}
