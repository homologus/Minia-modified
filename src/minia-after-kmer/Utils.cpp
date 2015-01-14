#include "Utils.h"
#include "Bank.h"

// some globals that don't really belong anywhere
int nks; // min abundance
uint32_t max_couv = 2147483646; // note: uint_abundance_t is 32 bits in SortingCount.cpp
 struct timeval tim;

const char *solid_kmers_file = (char *)"solid_kmers_binary"; 
const char *false_positive_kmers_file = (char *)"false_positive_kmers";
const char *assembly_file = (char *)"contigs.fa";
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
        sprintf(fileName,"%s.%s",prefix,suffix);
    else
        sprintf(fileName,"%s",suffix);
    return fileName;
}


int readlen;


int estimated_BL1;
uint64_t estimated_BL1_freesize;

float NBITS_PER_KMER = 11 ; // number of bits per kmer that optimizes bloo1 size


// below this line: unused kmer counting code

FILE * F_kmercpt_read;
FILE * F_kmercpt_write;





void Progress::init(uint64_t ntasks, const char * msg)
{
    gettimeofday(&timestamp, NULL);
    heure_debut = timestamp.tv_sec +(timestamp.tv_usec/1000000.0);
    
    fprintf(stderr,"| %-*s |\n",98,msg);
    
    todo= ntasks;
    done = 0;
    partial =0;
    for (int ii=0; ii<16;ii++) partial_threaded[ii]=0;
    for (int ii=0; ii<16;ii++) done_threaded[ii]=0;
    subdiv= 100;
    steps = (double)todo / (double)subdiv;
    
    if(!timer_mode)
    {
        fprintf(stderr,"[");fflush(stderr);
    }
}

void Progress::finish()
{
    set(todo);
    if(timer_mode)
        fprintf(stderr,"\n");
    else
        fprintf(stderr,"]\n");
    
    fflush(stderr);
    todo= 0;
    done = 0;
    partial =0;
    
}


void Progress::finish_threaded()// called by only one of the threads
{
    done = 0;
    double rem = 0;
    for (int ii=0; ii<16;ii++) done += (done_threaded[ii] ); 
    for (int ii=0; ii<16;ii++) partial += (partial_threaded[ii] ); 

    finish();

}

void Progress::inc(uint64_t ntasks_done, int tid)
{
    partial_threaded[tid] += ntasks_done;
    done_threaded[tid] += ntasks_done;
    while(partial_threaded[tid] >= steps)
    {
        if(timer_mode)
        {
            struct timeval timet;
            double now;
            gettimeofday(&timet, NULL);
            now = timet.tv_sec +(timet.tv_usec/1000000.0);
            uint64_t total_done  = 0;
            for (int ii=0; ii<16;ii++) total_done += (done_threaded[ii] );
            double elapsed = now - heure_debut;
            double speed = total_done / elapsed;
            double rem = (todo-total_done) / speed;
            if(total_done > todo) rem =0;
            int min_e  =  (int)(elapsed / 60) ;
            elapsed -= min_e*60;
            int min_r  =  (int)(rem / 60) ;
            rem -= min_r*60;
            
            fprintf(stderr,"%c%-5.3g  %%     elapsed: %6i min %-4.0f  sec      estimated remaining: %6i min %-4.0f  sec ",13,100*(double)total_done/todo,min_e,elapsed,min_r,rem);

        }
        else
        {
            fprintf(stderr,"-");fflush(stderr);
        }
        partial_threaded[tid] -= steps;

    }
    
}

void Progress::inc(uint64_t ntasks_done)
{
    done += ntasks_done;
    partial += ntasks_done;
    
    
    while(partial >= steps)
    {
        if(timer_mode)
        {
            gettimeofday(&timestamp, NULL);
            heure_actuelle = timestamp.tv_sec +(timestamp.tv_usec/1000000.0);
            double elapsed = heure_actuelle - heure_debut;
            double speed = done / elapsed;
            double rem = (todo-done) / speed;
            if(done>todo) rem=0;
            int min_e  = (int)(elapsed / 60) ;
            elapsed -= min_e*60;
            int min_r  = (int)(rem / 60) ;
            rem -= min_r*60;
            
            fprintf(stderr,"%c%-5.3g  %%     elapsed: %6i min %-4.0f  sec      estimated remaining: %6i min %-4.0f  sec ",13,100*(double)done/todo,min_e,elapsed,min_r,rem);
        }
        else
        {
            fprintf(stderr,"-");fflush(stderr);
        }
        partial -= steps;
    }
    
    
}


void Progress::set(uint64_t ntasks_done)
{
    if(ntasks_done > done)
        inc(ntasks_done-done);
}


