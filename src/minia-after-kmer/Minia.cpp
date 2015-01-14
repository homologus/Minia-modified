#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <vector> // for sorting_kmers
#include <sys/time.h>

#define NNKS 4 // default minimal abundance for solidity
#define MIN_CONTIG_SIZE 108

int max_memory; // the most memory one should alloc at any time, in MB

int order = 0; // deblooming order; 0 = debloom everything; 1 = don't debloom 1-node tips (experimental, untested, shouldn't work);// (made extern int in Traversal.h)
  
#include "Bank.h"
#include "Utils.h"
#include "SortingCount.h"
#include "Kmer.h"

int64_t genome_size;


int main(int argc, char *argv[])
{
    
    if(argc <  6)
    {
        fprintf (stderr,"usage:\n");
        fprintf (stderr," %s fasta_file kmer_size min_abundance estimated_genome_size prefix\n",argv[0]);
        fprintf (stderr,"hints:\n min_abundance ~ 3\n estimated_genome_size is in bp, does not need to be accurate, only controls memory usage\n prefix is any name you want the results to start with\n");

        return 1;
    }

    bool FOUR_BLOOM_VERSION = true;

     // shortcuts to go directly to assembly using serialized bloom and serialized hash
    int START_FROM_SOLID_KMERS=0; // if = 0, construct the fasta file of solid kmers, if = 1, start directly from that file 
    int LOAD_FALSE_POSITIVE_KMERS=0; // if = 0, construct the fasta file of false positive kmers (debloom), if = 1, load that file into the hashtable
    int NO_FALSE_POSITIVES_AT_ALL=0; // if = 0, normal behavior, if = 1, don't load false positives (will be a probabilistic de bruijn graph)
    int max_disk_space = 0;// let dsk decide
    for (int n_a = 6; n_a < argc ; n_a++)
    {
        if (strcmp(argv[n_a],"--original") == 0)
    	    FOUR_BLOOM_VERSION = false;

        if (strcmp(argv[n_a],"--dont-count")==0)
            START_FROM_SOLID_KMERS = 1;

        if (strcmp(argv[n_a],"--dont-debloom")==0)
            LOAD_FALSE_POSITIVE_KMERS = 1;

        if (strcmp(argv[n_a],"--just-assemble")==0)
        {
            START_FROM_SOLID_KMERS = 1;
            LOAD_FALSE_POSITIVE_KMERS = 1;
        }

        if (strcmp(argv[n_a],"--titus-mode")==0)
            NO_FALSE_POSITIVES_AT_ALL = 1;
        
        
        if (strcmp(argv[n_a],"-d")==0)
            max_disk_space = atoi(argv[n_a+1]);
        
        
        if (strcmp(argv[n_a],"-maxc")==0)
	    max_couv = atoi(argv[n_a+1]);
        
        if (strcmp(argv[n_a],"--le-changement")==0)
            {printf("c'est maintenant!\n");exit(0);}
    }


    // kmer size
    sizeKmer=27; // let's make it even for now, because i havnt thought of how to handle palindromes (dont want to stop on them)
    if(argc >=  3)
    {
        sizeKmer = atoi(argv[2]);
        if (sizeKmer%2==0)
        {
            sizeKmer-=1;
            printf("Need odd kmer size to avoid palindromes. I've set kmer size to %d.\n",sizeKmer);
        }
        if (sizeKmer>((int)sizeof(kmer_type)*4))
        {
            printf("Max kmer size on this compiled version is %d\n",sizeof(kmer_type)*4);
            exit(1);
        }
    }

    kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
    double lg2 = log(2);
   
    if (sizeKmer > 128)
    {
        FOUR_BLOOM_VERSION = false;
        printf("Reverted to single Bloom filter implementation for k>128\n");
    }

    // solidity 
    nks =NNKS;
    if(argc >=  4)
    {
        nks = atoi(argv[3]);
        if (nks==0) nks=1; // min abundance can't be 0
    }


   if(argc >=  5)
    {
       genome_size  = atoll(argv[4]);
      // int estimated_bloom_size = max( (int)ceilf(log2f(genome_size * NBITS_PER_KMER )), 1);
        uint64_t estimated_bloom_size = (uint64_t) genome_size * NBITS_PER_KMER;

       uint64_t estimated_nb_FP =  (uint64_t)(genome_size * 4 * powf(0.6,11)); // just indicative
    
       //max_memory = max( (1LL << estimated_bloom_size)/8LL /1024LL/1024LL, 1LL );
        max_memory =  max((int64_t) estimated_bloom_size/8LL /1024LL/1024LL,1LL);

      printf("estimated values: nbits Bloom %lli, nb FP %lld, max memory %i MB\n",estimated_bloom_size,estimated_nb_FP,max_memory);

    }

    // output prefix
    if(argc >=  6)
    {
        strcpy(prefix,argv[5]);
    }


    STARTWALL(0);

    Bank *Reads = new Bank(argv[1]);
    
    // counter kmers, write solid kmers to disk
    if (!START_FROM_SOLID_KMERS)
    {
        int verbose = 0;
        bool write_count = false;

        sorting_count(Reads,prefix,max_memory,max_disk_space,write_count,verbose);
    }

    return 0;
}


