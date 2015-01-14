#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <string>
#include <unistd.h>
#include <sys/types.h>
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <sys/time.h>
#include <math.h>

#include "minia/Bank.h"
#include "minia/Kmer.h"
#include "minia/Utils.h"
#include "minia/SortingCount.h" // for uint_abundance_t

using namespace std;

// some global variables..
uint64_t max_memory; // the most memory one should alloc at any time, in MB
int coverage_threshold;
int max_tip_length;

char *solid_kmers_with_count_filename;
char out_filename[1024];
OAHash *hash;
long nb_nodes_first_pass = 0, nb_nodes_second_pass = 0;
long nb_tips_removed = 0, nb_nodes_kept = 0;
bool hash_empty;
long nb_nodes_since_last_batch = 0;
FILE * file_new_nodes;


struct __attribute__((__packed__)) kmer_count_t // avoid any padding because i'll memcpy to it
{
    kmer_type kmer;
    uint_abundance_t abundance;
} kmer_count;

void stream_kmer_counts_and_populate_hashtable()
{
    BinaryBank *SolidKmersWithCount = new BinaryBank(solid_kmers_with_count_filename,sizeof(kmer_count_t),false);
    uint64_t header;
    SolidKmersWithCount->read(&header, 8);
    kmer_count_t kmer_count;

    while (SolidKmersWithCount->read_element_buffered(&kmer_count))
    {
        if (hash->has_key(kmer_count.kmer))
        {
            hash->insert(kmer_count.kmer, kmer_count.abundance);
        }
    }
    SolidKmersWithCount->close();
} 

void stream_nodes_second_pass(Bank *Nodes, long nb_nodes_to_read)
{
    long nodes_read_so_far = 0;
    bool eof = false;
    while(1)
    {
        char * seq;
        int seqlen;

        if(! Nodes->get_next_seq(&seq,&seqlen))
        {
            eof = true;
            break;
        }

        if (seqlen == 0) continue;
        bool is_tip = false;
        float avgcov = 0;
        int nbkmers = seqlen - sizeKmer + 1;
        bool has_left_neighbor = false, has_right_neighbor = false;

        bool is_short_seq = (seqlen <= max_tip_length);
        if (is_short_seq)
        {
            is_tip = true;
            kmer_type kmer, graine, graine_revcomp;

            for (int j=0; j<nbkmers; j++)
            {
                kmer = extractKmerFromRead(seq,j,&graine,&graine_revcomp, true, sizeKmer, kmerMask);

                int abundance_int;
                bool ret = hash->get(kmer,&abundance_int);
                if (!ret)
                {
                    printf("Unexpected error. %ld/%ld kmer %s from a sequence of length %d\n",nodes_read_so_far,nb_nodes_to_read,print_kmer(kmer),seqlen);
                    exit(1);
                }

                // was another criteria that i'm not using anymore.. it's not a tip if just one k-mer is above coverage_threshold
                /*if (abundance_int >= coverage_threshold)
                    is_tip = false;*/

                // test presence of left neighbors
                if (j == 0)
                {
                    for(int nt=0; nt<4; nt++) 
                    {
                        int current_strand = 1; 
                        kmer_type neighbor_kmer = next_kmer(kmer,nt,&current_strand);
                        int neighbor_abundance_int;
                        hash->get(neighbor_kmer,&neighbor_abundance_int);
                        has_left_neighbor |= (neighbor_abundance_int != -1);
                    }
                }

                // test presence of right neighbors
                if (j == nbkmers-1)
                {
                    for(int nt=0; nt<4; nt++) 
                    {
                        int current_strand = 0;
                        kmer_type neighbor_kmer = next_kmer(kmer,nt,&current_strand);
                        int neighbor_abundance_int;
                        hash->get(neighbor_kmer,&neighbor_abundance_int);
                        has_right_neighbor |= (neighbor_abundance_int != -1);
                    }
                }

                avgcov += abundance_int;
            }
            avgcov /= nbkmers;
        }

        if (avgcov > coverage_threshold) 
            is_tip = false;

        if (has_left_neighbor && has_right_neighbor)
            is_tip = false;

        // write seq if it's not a tip
        if (! is_tip)
        {
            if (is_short_seq)
                fprintf(file_new_nodes,">seq%ld__avgcov__%0.2f%s%s\n",nb_nodes_second_pass, avgcov, has_left_neighbor?"__has_left_neighbor":"", has_right_neighbor?"__has_right_neighbor":"");
            else
                fprintf(file_new_nodes,">seq%ld\n",nb_nodes_second_pass);
            fprintf(file_new_nodes,"%s\n",seq);
            nb_nodes_kept++;
        }
        else
            nb_tips_removed++;  

        nodes_read_so_far++;
        nb_nodes_second_pass++;
        if (nodes_read_so_far == nb_nodes_to_read)
            break;
    }

    if (!eof)
        printf("Processing node %d...\n",nb_nodes_second_pass); 
}

void process_batch(Bank *Nodes)
{
    stream_kmer_counts_and_populate_hashtable();

    //printf("processing batch (%d reads)\n",nb_nodes_first_pass - nb_nodes_since_last_batch);
    Nodes->load_position();
    stream_nodes_second_pass(Nodes, nb_nodes_first_pass - nb_nodes_since_last_batch);
    Nodes->save_position();

    nb_nodes_since_last_batch = nb_nodes_first_pass;

    delete hash;
    hash = new OAHash(max_memory);
    hash_empty = true;

    // at this point we're ready to continue processing the next batch of nodes
}


void stream_nodes_main_pass(Bank *Nodes)
{
    hash = new OAHash(max_memory);
    Nodes->save_position();
    file_new_nodes = fopen(out_filename,"w");
    hash_empty = true;

    while(1)
    {
        if (hash->load_factor() > 0.7)
            process_batch(Nodes);

        char * seq;
        int seqlen;

        if(! Nodes->get_next_seq(&seq,&seqlen)) break;
        if (seqlen == 0) continue;

        if (seqlen <= max_tip_length)
        {
            int nbkmers = seqlen - sizeKmer + 1;
            kmer_type kmer, graine, graine_revcomp;

            for (int j=0; j<nbkmers; j++)
            {
                kmer = extractKmerFromRead(seq,j,&graine,&graine_revcomp, true, sizeKmer, kmerMask);

                hash->insert(kmer, -1);
                hash_empty = false;

                // mark possible left neighbors for later examination
                if (j == 0)
                {
                    for(int nt=0; nt<4; nt++) 
                    {
                        int current_strand = 1;
                        kmer_type neighbor_kmer = next_kmer(kmer,nt,&current_strand);
                        hash->insert(neighbor_kmer, -1);
                    }
                }

                // mark possible right neighbors for later examination
                if (j == nbkmers-1)
                {
                    for(int nt=0; nt<4; nt++) 
                    {
                        int current_strand = 0;
                        kmer_type neighbor_kmer = next_kmer(kmer,nt,&current_strand);
                        hash->insert(neighbor_kmer, -1);
                    }
                }
            }
        }

        nb_nodes_first_pass++;
    }


    if (! hash_empty)
        process_batch(Nodes); // do the last batch

    printf("Finished. %ld tips removed, %ld non-tip nodes written.\n", nb_tips_removed, nb_nodes_kept);
    fclose(file_new_nodes);
}

int main(int argc, char *argv[])
{
    
    if(argc <  3)
    {
        fprintf (stderr,"parameters:\n");
        fprintf (stderr," nodes_file dsk_solid_kmers_file coverage_threshold [max_memory_in_MB]\n");
        return 0;
    }

    // first arg: input nodes file
    Bank *Nodes= new Bank(argv[1]);

    // second arg: kmer counts from DSK
    solid_kmers_with_count_filename = argv[2];

    // third arg: threshold for coverage
    coverage_threshold = atoi(argv[3]);
    if (coverage_threshold < 1)
    {
        printf("Coverage threshold must be > 0.\n");
        exit(1);
    }

    // output nodes file
    sprintf(out_filename,"%s.tip_filtered_nodes",argv[1]);

    // estimate two things from the solid k-mers file: k and the number of distinct kmers
    BinaryBank *SolidKmersWithCount = new BinaryBank(solid_kmers_with_count_filename,sizeof(kmer_count_t),false);
    uint32_t kmer_nbits;
    SolidKmersWithCount->read(&kmer_nbits, 4);
    uint32_t k;
    SolidKmersWithCount->read(&k, 4);
    SolidKmersWithCount->close();

    // set k-mer size now that we know it (and mask: Minia stuff)
    sizeKmer = k;
    max_tip_length = sizeKmer*2;
    kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;

    if (kmer_nbits/2 > (sizeof(kmer_type)*4))
    {
       printf("The solid k-mers file was created with `make k=%d` but this software was compiled with `make k=%d`\n", kmer_nbits/2, sizeof(kmer_type)*4);
       exit(1);
    }

    // estimate number of kmers in the DSK file (will just be off by 1 normally because of the header)
    long rough_nb_kmers = SolidKmersWithCount->nb_elements();
    printf("Number of distinct k-mers: %ld M\n",rough_nb_kmers/1024/1024);

    // set max memory based on fourth argument or 2 bits per kmer
    if (argc > 4)
        max_memory = max( atoll(argv[4]), 1LL) * 1024LL*1024LL; // max(specified memory, 1 MB)
    else
        max_memory = max( (rough_nb_kmers * 2)/8LL, 1024LL*1024LL); // max(2*(genome size), 1 MB)

    printf("Maximum memory: %d MB\n",(int)(max_memory/1024LL/1024LL));

    // now actually start the computation (timed)

    STARTWALL(0);

    printf("Removing tips from \"%s\", writing results to \"%s\"\n",argv[1],out_filename);
    printf("Criteria for removal: node length <= 2k (= %d bp) and mean k-mer coverage <= %d\n",max_tip_length,coverage_threshold);

    stream_nodes_main_pass(Nodes);
    
    STOPWALL(0,"Total");

    return 0;
}


