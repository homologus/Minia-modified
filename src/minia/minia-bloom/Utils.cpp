#include "Utils.h"
#include "Bank.h"

// some globals that don't really belong anywhere
int nks; // min abundance
uint32_t max_couv = 2147483646; // note: uint_abundance_t is 32 bits in SortingCount.cpp
 struct timeval tim;

const char *solid_kmers_file = (char *)"solid_kmers_binary"; 
const char *false_positive_kmers_file = (char *)"false_positive_kmers";
const char *bloom_file = (char *)"bloom_data";
const char *assembly_file = (char *)"contigs.fa";
const char *branching_kmers_file = (char *)"branching_kmers"; // (only useful for multiple assemblies with same bloom&debloom structure (ie debugging))
const char *binary_read_file = (char *)"reads_binary";
const char *histo_file_name = (char *)"histo";
const char *breakpoints_file_name = (char *)"breakpoints";

const char *assoc_kmer_file = (char *)"paired_kmer";


// prefix-based output files naming 
char prefix[1024];
char fileName[1024];
char *return_file_name(const char *suffix)
{
    if (strlen(prefix)>0)
        sprintf(fileName,"a.%s",suffix);
    else
        sprintf(fileName,"a.%s",suffix);
    return fileName;
}


int readlen;

template <typename T, typename U>
// T can be Bloom, BloomCpt, BloomCpt3 or LinearCounter (just needs to support add(kmer_type) and possibly contains(kmer_type))
// U can be BloomCpt or BloomCpt3
void bloom_pass_reads(Bank *Sequences, T *bloom_to_insert, U *bloom_counter, char *stderr_message)
{
    int64_t NbRead = 0;
    int64_t NbInsertedKmers = 0;
    Sequences->rewind_all();
    char * rseq;
    long i;
    kmer_type kmer, graine, graine_revcomp;


    while (Sequences->get_next_seq(&rseq,&readlen))
    {
      for (i=0; i<readlen-sizeKmer+1; i++)
        {
            kmer = extractKmerFromRead(rseq,i,&graine,&graine_revcomp);

            if (bloom_counter != NULL)
            {
                // discard kmers which are not solid
                if( ! bloom_counter->contains_n_occ(kmer,nks)) continue;
            }

            bloom_to_insert->add(kmer);
            NbInsertedKmers++;

        }
        NbRead++;
        if ((NbRead%10000)==0) fprintf (stderr,stderr_message,13,NbRead);
    }
    fprintf (stderr,"\nInserted %lld %s kmers in the bloom structure.\n",(long long)NbInsertedKmers,"(redundant)");

}


template <typename T> // T can be Bloom, BloomCpt or BloomCpt3
void bloom_pass_reads_binary(T *bloom_to_insert, BloomCpt *bloom_counter, char *stderr_message)
{
  fprintf(stderr,"binary pass \n");
  int64_t NbRead = 0;
  int64_t NbInsertedKmers = 0;
  kmer_type kmer;
  
  // read solid kmers from disk
  BinaryBank * SolidKmers = new BinaryBank(return_file_name(solid_kmers_file),sizeof(kmer),0);

  while(SolidKmers->read_element(&kmer))
    {
      // printf("kmer %lld\n",kmer);
      bloom_to_insert->add(kmer);
      NbInsertedKmers++;
      NbRead++;
      if ((NbRead%10000)==0) fprintf (stderr,stderr_message,13,(long long)NbRead);
    }
  fprintf (stderr,"\nInserted %lld %s kmers in the bloom structure.\n",(long long)NbInsertedKmers,"solid");
  SolidKmers->close();
  
}

int estimated_BL1;
uint64_t estimated_BL1_freesize;

float NBITS_PER_KMER = 11 ; // number of bits per kmer that optimizes bloo1 size



// loading bloom from disk
template <typename T> //bloocpt or bloocpt3
Bloom *bloom_create_bloo1(T *bloom_counter, bool from_dump)
{

    BinaryBank * SolidKmers ;
    
    if(from_dump && nsolids) // from dump and known number of solid kmers 
    {
        //nsolids is sotred in a config file
        //number of solid kmers cannot be computed precisely from bloom file, imprecision of 0-7
        estimated_BL1 = max( (int)ceilf(log2f(nsolids*NBITS_PER_KMER)), 1);
        estimated_BL1_freesize =  (uint64_t)(nsolids*NBITS_PER_KMER);
    }
    else
    {
        // get true number of solid kmers, in order to precisely size the bloom filter
        SolidKmers = new BinaryBank("a.solid_kmers_binary",sizeof(kmer_type),0);
        estimated_BL1 = max( (int)ceilf(log2f(SolidKmers->nb_elements()*NBITS_PER_KMER)), 1);
        estimated_BL1_freesize =  (uint64_t)(SolidKmers->nb_elements()*NBITS_PER_KMER);
        printf("nelem %lli nbits %g \n",SolidKmers->nb_elements(),NBITS_PER_KMER);
    }
    
    //printf("Allocating %0.1f MB of memory for the main Bloom structure (%g bits/kmer)\n",(1LL<<estimated_BL1)/1024.0/1024.0/8.0,NBITS_PER_KMER);
    printf("freesize %lli estimated_BL1_freesize  %0.1f MB of memory for the main Bloom structure (%g bits/kmer)\n",(long long)estimated_BL1_freesize,(estimated_BL1_freesize)/1024.0/1024.0/8.0,NBITS_PER_KMER);
    
    Bloom *bloo1;
#if CUSTOMSIZE
    bloo1 = new Bloom((uint64_t)estimated_BL1_freesize);
#else
    bloo1 = new Bloom(estimated_BL1);
#endif

    bloo1->set_number_of_hash_func((int)floorf(0.7*NBITS_PER_KMER));

        bloom_pass_reads_binary(bloo1, bloom_counter, (char*)"%cInsert solid Kmers in Bloom %lld"); // use the method reading SolidKmers binary file, was useful when varying Bloom size (!= dumped size)
        //bloo1->dump(return_file_name(bloom_file)); // create bloom dump
        SolidKmers->close();

    return bloo1;    
}

// wrapper for default behavior: don't load from dump
template <typename T> //bloocpt or bloocpt3
Bloom *bloom_create_bloo1(T *bloom_counter)
{
    return bloom_create_bloo1(bloom_counter, false);
}

template Bloom *bloom_create_bloo1<BloomCpt>(BloomCpt *bloom_counter); // trick to avoid linker errors: http://www.parashift.com/c++-faq-lite/templates.html#faq-35.13
template Bloom *bloom_create_bloo1<BloomCpt>(BloomCpt *bloom_counter, bool from_dump); 


