## k-mer frequency analysis

### The basic use of k-mer

- 1) Genome estimation
- 2) Basic units for assembling
- 3) Exploring sequence composition

### kmer_count.pl
```
 Script for kmer counting.
	perl kmer_count.pl seq.fna [options]

	--AsOne            -a  if multi-sequences in a fasta file, count as one

	--ksize            -k  k-mer size, default is 4

	--canonical_kmer   -c  canonical k-mer type
```


### Test based on TNF (Tetranucleotide frequency)

<div align=center><img src="https://github.com/wangpeng407/similarity_based_on_k-mer/blob/master/otu_kmer_test/pca_tnf.png"/></div>
