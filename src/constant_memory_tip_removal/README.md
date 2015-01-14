This software removes tips and short isolated vertices from a de Bruijn graph that resides in a file, in constant memory

input:
 - set of nodes (fasta)
 - kmer counts (output of DSK for the reads)
 - coverage threshold

output: 
- filtered set of nodes (fasta)

definitions:
- a tip is a graph node such that:
  * it has no in-neighbors xor no out-neighbors
  * mean k-mer coverage < coverage_threshold (this definition is my own, Velvet for example only examines the coverage of the edge that leads to the tip)
  * node length < max_tip_length
- an isolated vertex is a node having no in-neighbor and no out-neighbor and has length < max_tip_length
- max_tip_length is set to 2k (similar to Velvet)
- coverage_threshold is a user-defined parameter

For convenience, the default memory usage is set to be (2*number of distinct k-mers) bits

known bugs:
- if a node is larger than (genome_size/4) nucleotides then this software will use more than (2*genome_size) memory

license: CECILL (uses parts of Minia)

author: Rayan Chikhi

September 2013
