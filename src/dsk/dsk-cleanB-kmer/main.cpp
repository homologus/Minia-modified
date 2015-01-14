#include "SortingCount.h"
#include "inttypes.h"
#include <sys/resource.h> // for getrlimit()

#define SINGLE_BAR 1

bool clear_cache = false; // clear file cache from memory (for timing only)

bool hybrid_mode = false;
bool use_hashing = true; // use hashing instead of sorting (better control of memory)
float load_factor = 0.7;
bool use_compressed_reads = true ; // true; // write compressed read file

bool output_histo;


// main k-mer counting function, shared between minia and dsk
// verbose == 0 : stderr progress bar
// verbose >= 1 : print basic status
// verbose >= 2 : print extra partition information
// write_count == True: include kmer count in results file, in that form:
//           - save kmer count for each kmer in the resulting binary file
//           - the very first four bytes of the result file are the kmer length
void sorting_count(Bank *Sequences, char *prefix, int max_memory, int max_disk_space, bool write_count, int verbose)
{

    // create a temp dir from the prefix
    char temp_dir[1024];
    sprintf(temp_dir,"%s_temp",prefix);

    // clear the temp folder (needs to be done before estimating disk space)
    DIR*            dp;
    struct dirent*  ep;
    char            p_buf[512] = {0};
    dp = opendir(temp_dir);
    while ( (dp != NULL) && ((ep = readdir(dp)) != NULL)) {
        sprintf(p_buf, "%s/%s", temp_dir, ep->d_name);
        remove(p_buf);
    }
    if(dp != NULL)
        closedir(dp);

    if (max_disk_space == 0)
    {
        // default max disk space
        struct statvfs buffer ;
        char current_path[1000];
        getcwd(current_path,sizeof(current_path));
        // int ret =
        statvfs(current_path, &buffer);
        int available = (int)(((double)buffer.f_bavail * (double)buffer.f_bsize) / 1024 / 1024);
        printf("Available disk space in %s: %d MB\n",current_path,available); // not working in osx (is that a TODO then?)
        max_disk_space = min(available/2, (int)(( (double)Sequences->filesizes ) / 1024 / 1024));
    } 
    if (max_disk_space == 0) // still 0?
        max_disk_space = 10000; // = default for osx

    // estimate number of iterations
    uint64_t volume = Sequences->estimate_kmers_volume(sizeKmer);
    uint32_t nb_passes = ( volume / max_disk_space ) + 1;

    
    int nb_threads=1;
    
    
    // temp bugfix: don't use compressed reads for long reads
    if (Sequences->estimate_max_readlen() > 1000000)
        use_compressed_reads = false;
    
    
    uint64_t volume_per_pass;
    uint32_t nb_partitions;


    // loop to lower the number of partitions below the maximum number of simulatenously open files
    do
    {
        volume_per_pass = volume / nb_passes;
        nb_partitions = ( volume_per_pass / max_memory ) + 1;

        // if partitions are hashed instead of sorted, adjust for load factor
        // (as in the worst case, all kmers in the partition are distinct and partition may be slightly bigger due to hash-repartition)
        if (use_hashing)
        {
            nb_partitions = (uint32_t) ceil((float) nb_partitions / load_factor);
            nb_partitions = ((nb_partitions * OAHash::size_entry()) + sizeof(key_type)-1) / sizeof(key_type); // also adjust for hash overhead
            //printf("Updated number of partitions for hash-based k-mer counting: %d\n",nb_partitions);
        }

        // round nb_partitions to mulitple of nthreads, for better perf
        //  nb_partitions = ((nb_partitions + nb_threads - 1) / nb_threads) * nb_threads;
        
        // get max number of open files
        struct rlimit lim;
        int max_open_files = 1000;
        int err = getrlimit(RLIMIT_NOFILE, &lim);
        if (err == 0)
            max_open_files = lim.rlim_cur / 2;
        if (nb_partitions >= max_open_files)
            nb_passes++;
        else
            break;
    }
    while (1);

 // volume / (sizeof(kmer_type)*4)   is approx size of read file stored in binary, read nb_passes -1 times
    uint64_t total_IO =   volume * 2LL * 1024LL*1024LL   ;// in bytes  +   nb_passes * ( volume / (sizeof(kmer_type)*4) )    ; // in bytes
    uint64_t temp_IO = 0;
    //if (nb_passes==1) use_compressed_reads=false;
    BinaryBankConcurrent * redundant_partitions_file[nb_partitions]; 
    char redundant_filename[nb_partitions][256];
    kmer_type kmer;
    int max_read_length = KMERSBUFFER_MAX_READLEN;
    kmer_type * kmer_table_seq = (kmer_type * ) malloc(sizeof(kmer_type)*max_read_length); ;


    BinaryReads *  binread = NULL;
    if(use_compressed_reads)
        binread = new BinaryReads(return_file_name(binary_read_file),true);

    fprintf(stderr,"Sequentially counting ~%llu MB of kmers with %d partition(s) and %d passes using %d thread(s), ~%d MB of memory and ~%d MB of disk space\n", (unsigned long long)volume, nb_partitions,nb_passes, nb_threads, max_memory * nb_threads, max_disk_space);

    STARTWALL(count);

    mkdir(temp_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    BinaryBankConcurrent * SolidKmers = new BinaryBankConcurrent(return_file_name(solid_kmers_file),sizeof(kmer),true,nb_threads);

    if (write_count)
    {
        // write k-mer nbits as the first 4 bytes; and actual k-mer size as the next 4 bits
        uint32_t kmer_nbits = sizeof(kmer) * 8;
        SolidKmers->write_buffered(&kmer_nbits, 4,0);
        SolidKmers->write_buffered(&sizeKmer, 4,0);
        SolidKmers->flush(0);
    }

    int64_t estimated_NbReads = Sequences->estimate_nb_reads();
    char * rseq;
    int readlen;
    int64_t NbSolid = 0;
    int64_t * NbSolid_omp = (int64_t  *) calloc(nb_threads,sizeof(int64_t));
    //long total_kmers_per_partition[nb_partitions]; //guillaume probably commented it because updating this variable would require synchronization
    long distinct_kmers_per_partition[nb_partitions];
    uint64_t  * histo_count = (uint64_t  *) calloc(10001,sizeof(uint64_t));


    //start by the conversion of the file to binary format

    if(use_compressed_reads)
    {
        char * pt_begin;
        int idx =0 ;
        int64_t NbRead = 0;
        Progress progress_conversion;
       // progress_conversion.timer_mode=1; // to switch to timer mode (show elapsed and estimated remaining time)
        progress_conversion.init(estimated_NbReads,"First step: Converting input file into Binary format");
        
        Sequences->rewind_all();
        while(1)
        {
            if(! Sequences->get_next_seq(&rseq,&readlen)) break; // read  original fasta file
            if(readlen > max_read_length) // realloc kmer_table_seq if needed
            {
                max_read_length = 2*readlen;
                kmer_table_seq = (kmer_type * ) realloc(kmer_table_seq,sizeof(kmer_type)*max_read_length);
            }
            
            pt_begin = rseq;
            
            //should be ok
            while (pt_begin < (rseq+ readlen))
            {
                idx=0; // start a new read

                //skips NN
                while (*pt_begin =='N' && pt_begin < (rseq+ readlen))
                {
                    pt_begin ++;
                }
                // goes to next N or end of seq
                while ( (pt_begin[idx] !='N') &&  ((pt_begin +idx) < (rseq+ readlen))  )
                {
                    idx++;
                }
                
                //we have a seq beginning at  pt_begin of size idx  ,without any N, will be treated as a read:
                binread->write_read(pt_begin,idx);
                pt_begin += idx;
            }
            
            // binread->write_read(rseq,readlen);
            
            
            NbRead++;
            if ((NbRead%10000)==0)
            {
                progress_conversion.inc(10000);
            }
        }
        progress_conversion.finish();
        binread->close();

    }
    ///fin conversion
    if (clear_cache)
    {
#ifdef OSX
        system("purge");
#else
        system("echo 3 > /proc/sys/vm/drop_caches");
#endif
    }
    
    
    
#if SINGLE_BAR
    Progress progress;
    char message[1000];
    sprintf(message,"Counting kmers");
    progress.timer_mode=1;
    if (verbose == 0 )
        progress.init(total_IO,message);
#endif
    
    
    // how many times we will traverse the whole reads file (has an influence on temp disk space)
    for (uint32_t current_pass = 0; current_pass < nb_passes; current_pass ++)
    {
        
        if(use_compressed_reads ) //open binary reads for reading
            binread->open(false);
        
        STARTWALL(debpass);
        STARTWALL(debw);

        for (uint32_t p=0;p<nb_partitions;p++)
        {
            sprintf(redundant_filename[p],"%s/partition%d.redundant_kmers",temp_dir,p);
            redundant_partitions_file[p] =  new BinaryBankConcurrent (redundant_filename[p],sizeof(kmer_type),true, nb_threads);
            distinct_kmers_per_partition[p]=0;
        }

        // partitioning redundant kmers
        
        Sequences->rewind_all();
#if !SINGLE_BAR
        Progress progress;
        progress.timer_mode=1; // to switch to timer mode (show elapsed and estimated remaining time)
        char message[1000];
        sprintf(message,"Pass %d/%d, Step 1: partitioning",current_pass+1,nb_passes);
        if (verbose == 0 )
            progress.init(estimated_NbReads,message);
#endif
     

        
        //current_pass> 0 &&
        {
            int64_t  nbkmers_written =0;
            int tid =0;
            int64_t NbRead = 0;
            int64_t nread =0;
            int64_t tempread =0;
            int nreads_in_buffer= 1000;
            KmersBuffer * kbuff =NULL;
            if(use_compressed_reads)
            {
                kbuff = new KmersBuffer (binread, 1000000,  nreads_in_buffer); //buffer size (in nb of kmers), seq per task // the buffer is per thread
                kbuff->binary_read_file = binread->binary_read_file;
            }

            kmer_type * kmer_table ;
            while(1)
            {

                //read the fasta file
                if(use_compressed_reads) // && current_pass>0
                {
                    nread = kbuff->readkmers();
                    if(! nread) break;
                    NbRead+= nread;
                    tempread+= nread;
                }
                else
                {
                    if(! Sequences->get_next_seq(&rseq,&readlen)) break; // read  original fasta file
                    if(readlen > max_read_length) // realloc kmer_table_seq if needed
                    {
                        max_read_length = 2*readlen;
                        kmer_table_seq = (kmer_type * ) realloc(kmer_table_seq,sizeof(kmer_type)*max_read_length);
                    }

                }

//                if(use_compressed_reads ) //write compressed read file at first pass //&& current_pass==0
//                    binread->write_read(rseq,readlen);

                int i;
                int nbkmers =readlen-sizeKmer+1;

                if( use_compressed_reads) //current_pass >0 &&
                {
                    nbkmers = kbuff->nkmers;
                    kmer_table = kbuff->kmers_buffer;
                   // printf("nb kmers read  %lli \n",nbkmers);
                 //   NbRead+= nreads_in_buffer;
                } 
                else //old fashion
                {
                    compute_kmer_table_from_one_seq(readlen,rseq,kmer_table_seq);
                    nbkmers =readlen-sizeKmer+1;
                    kmer_table = kmer_table_seq;
                    NbRead++;
                }

                nbkmers_written= 0;
                //compute the kmers stored in the buffer kmer_table
                for (i=0; i<nbkmers; i++)
                {
                    kmer_type lkmer;

                    // kmer = extractKmerFromRead(rseq,i,&graine,&graine_revcomp);

                    lkmer = kmer_table[i];

                    // some hashing to uniformize repartition
                    kmer_type kmer_hash = lkmer ^ (lkmer >> 14);
                    kmer_hash = (~kmer_hash) + (kmer_hash << 18); 
                    kmer_hash = kmer_hash ^ (kmer_hash >> 31);
                    kmer_hash = kmer_hash * 21; 
                    kmer_hash = kmer_hash ^ (kmer_hash >> 11);
                    kmer_hash = kmer_hash + (kmer_hash << 6);
                    kmer_hash = kmer_hash ^ (kmer_hash >> 22);

                    // check if this kmer should be included in the current pass
                    if ((kmer_hash % nb_passes  ) != current_pass) 
                        continue;

                    kmer_type reduced_kmer = kmer_hash / nb_passes;

                    int p;// compute in which partition this kmer falls into

#ifdef _ttmath
                    (reduced_kmer % nb_partitions).ToInt(p);
#else
                    p = reduced_kmer % nb_partitions;
#endif

                    nbkmers_written++;

                    redundant_partitions_file[p]->write_element_buffered(&lkmer,tid); // save this kmer to the right partition file
                    // total_kmers_per_partition[p]++; // guillaume probably commented it because updating this variable would require synchronization

                }
                //NbRead++;
#if SINGLE_BAR
                if(verbose==0)
                {
                if (nb_threads == 1)
                    progress.inc(nbkmers_written * sizeof(kmer_type));
                else
                    progress.inc(nbkmers_written * sizeof(kmer_type),tid);
                }
#endif
             //   if ((NbRead%10000)==0)
                if(tempread> 10000)
                {
                    tempread -= 10000;
                    if (verbose)
                        fprintf (stderr,"%cPass %d/%d, loop through reads to separate (redundant) kmers into partitions, processed %lluM reads out of %lluM",13,current_pass+1,nb_passes,(unsigned long long)(NbRead/1000/1000),(unsigned long long)(estimated_NbReads/1000/1000));
#if !SINGLE_BAR
                    else
                        if (nb_threads == 1)
                            progress.set(NbRead);
                        else
                            progress.inc(10000,tid);
#endif
                }
            } //end while

            if(use_compressed_reads)
                delete kbuff;
        } 


        
#if !SINGLE_BAR
        if (verbose == 0)
        {
            if (nb_threads == 1)
             progress.finish();
            else
              progress.finish_threaded();  // here only one thread
            
            sprintf(message,"Pass %d/%d, Step 2: computing kmer count per partition",current_pass+1,nb_passes);
            progress.init(nb_partitions+1,message);
        }
#endif
        
        if (verbose)fprintf(stderr,"\n");

        if (verbose >= 2)
        {
            STOPWALL(debw,"Writing redundant kmers");
        }
        STARTWALL(debtri);

        // close partitions and open them for reading

            for (uint32_t p=0;p<nb_partitions;p++)
            {
                redundant_partitions_file[p]->close();
                redundant_partitions_file[p]->open(false);
            }



        // for better timing: clear the file cache, since the partitions may still be in memory, that's unfair to low mem machines
        if (clear_cache)
        {
#ifdef OSX
            system("purge");
#else
            system("echo 3 > /proc/sys/vm/drop_caches");
#endif
        }


        //quick and dirty parall with omp, testing
        //todo if we want omp and histo : separate histo_count tab per thread that needs to be merged at the end
        // TODO to guillaume: remove that todo above, because it is done, right?
        // load, sort each partition to output solid kmers
        for (int p=0;p<nb_partitions;p++)
        {
            kmer_type lkmer;
            
            bool use_hashing_for_this_partition = use_hashing;
            if(hybrid_mode)
            {
              //  printf("max mem %i MB  ,   parti size %i MB\n",max_memory,(redundant_partitions_file[p]->nb_elements()*sizeof(kmer_type))/1024LL/1024LL);
                if(   (redundant_partitions_file[p]->nb_elements()*sizeof(kmer_type)) <  (max_memory*1024LL*1024LL) )
                    use_hashing_for_this_partition = false;
                else
                    use_hashing_for_this_partition = true;
            }
            int tid =0;
            
            if (use_hashing_for_this_partition)
            {
                // hash partition and save to solid file
                OAHash hash(max_memory*1024LL*1024LL);
                uint64_t nkmers_read=0;
                
                while (redundant_partitions_file[p]->read_element_buffered (&lkmer))
                {
                    hash.increment(lkmer);
                    nkmers_read++;
#if SINGLE_BAR
                    if(verbose==0 && nkmers_read==10000)
                    {
                        if (nb_threads == 1)
                            progress.inc(nkmers_read*sizeof(kmer_type));
                        else
                            progress.inc(nkmers_read*sizeof(kmer_type),tid);
                        nkmers_read=0;
                    }
#endif
                }
                
                //single bar
                
                
                if (verbose >= 2)
                    printf("Pass %d/%d partition %d/%d hash load factor: %0.3f\n",current_pass+1,nb_passes,p+1,nb_partitions,hash.load_factor());
                
                hash.start_iterator();
                while (hash.next_iterator())
                {
                    uint_abundance_t abundance = hash.iterator->value;
                    if(output_histo)
                    {
                        uint_abundance_t saturated_abundance;
                        saturated_abundance = (abundance >= 10000) ? 10000 : abundance;
                        //printf("histo_count 0 1  2 %i %i %i \n",histo_count[0],histo_count[1],histo_count[2]);
                        
                        histo_count[saturated_abundance]++;
                    }
                    if (abundance >= nks && abundance <= max_couv)
                    {
                        SolidKmers->write_element_buffered(&(hash.iterator->key),tid);
                        
                        NbSolid_omp[tid]++;
                        if (write_count)
                            SolidKmers->write_buffered(&abundance, sizeof(abundance),tid, false);
                        
                    }
                    distinct_kmers_per_partition[p]++;
                }
            }
            
            else
            {
                // sort partition and save to solid file
                vector < kmer_type > kmers;
                uint64_t nkmers_read=0;
                
                
                
                while (redundant_partitions_file[p]->read_element_buffered (&lkmer))
                {
                    kmers.push_back (lkmer);
                    nkmers_read++;
#if SINGLE_BAR
                    if(verbose==0 && nkmers_read==10000)
                    {
                        if (nb_threads == 1)
                            progress.inc(nkmers_read*sizeof(kmer_type));
                        else
                            progress.inc(nkmers_read*sizeof(kmer_type),tid);
                        nkmers_read=0;
                    }
#endif
                }
                
                
                sort (kmers.begin (), kmers.end ());
                
                kmer_type previous_kmer = *(kmers.begin ());
                uint_abundance_t abundance = 0;
                for (vector < kmer_type >::iterator it = kmers.begin (); it != kmers.end ();
                     it++)
                {
                    kmer_type current_kmer = *it;
                    
                    if (current_kmer == previous_kmer)
                        abundance++;
                    else
                    {
                        if(output_histo)
                        {
                            uint_abundance_t saturated_abundance;
                            saturated_abundance = (abundance >= 10000) ? 10000 : abundance;
                            histo_count[saturated_abundance]++;
                            
                        }
                        if (abundance >= nks  && abundance <= max_couv)
                        {
                            NbSolid_omp[tid]++;
                            SolidKmers->write_element_buffered(&previous_kmer,tid);
                            
                            if (write_count)
                                SolidKmers->write_buffered(&abundance, sizeof(abundance),tid, false);
                        }
                        abundance = 1;
                        distinct_kmers_per_partition[p]++;
                    }
                    previous_kmer = current_kmer;
                }
                
                //last kmer
                distinct_kmers_per_partition[p]++;
                if(output_histo)
                {
                    uint_abundance_t saturated_abundance;
                    saturated_abundance = (abundance >= 10000) ? 10000 : abundance;
                    histo_count[saturated_abundance]++;
                    
                }
                if (abundance >= nks && abundance <= max_couv)
                {
                    NbSolid_omp[tid]++;
                    SolidKmers->write_element_buffered(&previous_kmer,tid);
                    
                    if (write_count)
                        SolidKmers->write_buffered(&abundance, sizeof(abundance),tid, false);
                    
                }
                
            }
            
            
            if (verbose >= 1)
                fprintf(stderr,"%cPass %d/%d, loaded and sorted partition %d/%d, found %lld solid kmers so far",13,current_pass+1,nb_passes,p+1,nb_partitions,(long long)(NbSolid_omp[tid]));
            
            if (verbose >= 2)
                printf("\nPass %d/%d partition %d/%d %ld distinct kmers\n",current_pass+1,nb_passes,p+1,nb_partitions,/*total_kmers_per_partition[p],*/distinct_kmers_per_partition[p]);
            
#if !SINGLE_BAR
            if (verbose == 0 && nb_threads==1)
                progress.inc(1);
            else if (verbose == 0 && nb_threads>1)
                progress.inc(1,tid);
#endif
            
            
            
            redundant_partitions_file[p]->close();
            remove(redundant_filename[p]);
            
        } // end for partitions

        
#if !SINGLE_BAR
        if (verbose == 0 && nb_threads == 1)
            progress.finish();
        else if (verbose == 0 && nb_threads > 1 )
            progress.finish_threaded();
#endif

        if (verbose) fprintf(stderr,"\n");

        if (verbose >= 2)
        {
            STOPWALL(debtri,"Reading and sorting partitions");
            STOPWALL(debpass,"Pass total");

        }
       
        if(use_compressed_reads)
            binread->close();
        
        //delete
            for (uint32_t p=0;p<nb_partitions;p++)
            {
                delete redundant_partitions_file[p] ;
            }
        
    }

    //single bar
#if SINGLE_BAR
    if (verbose == 0 && nb_threads == 1)
        progress.finish();
    else if (verbose == 0 && nb_threads > 1 )
        progress.finish_threaded();
#endif
    
    if(output_histo)
    {
        FILE * histo_file = fopen(return_file_name(histo_file_name),"w");
        for (int cc=1; cc<10001; cc++) {
            fprintf(histo_file,"%i\t%llu\n",cc,(unsigned long long)(histo_count[cc]));
        }
        fclose(histo_file);
    }
    free(histo_count);

    NbSolid = NbSolid_omp[0];

    SolidKmers->close();
    printf("\nSaved %lld solid kmers\n",(long long)NbSolid);
    rmdir(temp_dir);

    STOPWALL(count,"Counted kmers");
    fprintf(stderr,"\n------------------ Counted kmers and kept those with abundance >=%i,     \n",nks);
} 



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h> // for mkdir
#include <inttypes.h>
#include <stdint.h>
#include <algorithm> // for max/min
#include <vector> // for sorting_kmers
#include <sys/time.h>
#ifndef OSX
#include <sys/sysinfo.h> // to determine system memory
#endif
#include <sys/statvfs.h> // to determine available disk space
#include <dirent.h> // to clear the temp directory
#include <libgen.h> // for basename()
#include <string>

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

int max_memory; // the most memory dsk should alloc at any time, in MB
int max_disk_space; // the most disk space dsk should use at any time, in MB
extern bool output_histo;
#include "Bank.h"
#include "Utils.h"
#include "SortingCount.h"
#include "Kmer.h"

int main(int argc, char *argv[])
{
    if(argc <  3)
    {
        fprintf (stderr,"%s: [d]isk [s]treaming of [k]-mers (constant-memory k-mer counting)\n",argv[0]);
        fprintf (stderr,"usage:\n");
        fprintf (stderr," %s input_file kmer_size [-t min_abundance] [-m max_memory] [-d max_disk_space] [-o out_prefix] [-histo]\n",argv[0]);
        fprintf (stderr,"details:\n [-t min_abundance] filters out k-mers seen ( < min_abundance ) times, default: 1 (all kmers are returned)\n [-m max_memory] is in MB, default: min(total system memory / 2, 5 GB) \n [-d max_disk_space] is in MB, default: min(available disk space / 2, reads file size)\n [-o out_prefix] saves results in [out_prefix].solid_kmers. default out_prefix = basename(input_file)\n [-histo] outputs histogram of kmers abundance\n Input file can be fasta, fastq, gzipped or not, or a file containing a list of file names.\n");
#ifdef SVN_REV
fprintf(stderr,"Running dsk version %s\n",STR(SVN_REV));
#endif
        return 0;
    }

    // reads file
    Bank *Reads = new Bank(argv[1]);

    if (argv[2][0] == '-')
    {
        printf("please specify a k value\n");
        exit(1);
    }

    // kmer size
    sizeKmer = atoi(argv[2]);
    if (sizeKmer>(int)(sizeof(kmer_type)*4))
    {
        printf("Max kmer size on this compiled version is %lu\n",sizeof(kmer_type)*4);
        exit(1);
    }
    kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;

    // default solidity 
    nks = 1;

    // default max memory
    max_memory = 5*1024;
    #ifndef OSX
    struct sysinfo info;
    sysinfo(&info);
    int total_ram = (int)(((double)info.totalram*(double)info.mem_unit)/1024/1024);
    printf("Total RAM: %d MB\n",total_ram);
#else
    int total_ram = 128*1024;
#endif


    // default prefix is the reads file basename
    char *reads_path=strdup(argv[1]);
    string reads_name(basename(reads_path)); // posix basename() may alter reads_path
    free(reads_path);
    int lastindex = reads_name.find_last_of("."); 
    strcpy(prefix,reads_name.substr(0, lastindex).c_str()); 

    for (int n_a = 3; n_a < argc ; n_a++)
    {
        if (strcmp(argv[n_a],"-t")==0)
            nks = atoi(argv[n_a+1]);

        if (strcmp(argv[n_a],"-o")==0)
            strcpy(prefix,argv[n_a+1]);
    }

    int verbose = 0;

    max_disk_space = 0;

    output_histo =false;
    // parse the remaining arguments: these will override the default max memory / max disk
    for (int n_a = 3; n_a < argc ; n_a++)
    {
        if (strcmp(argv[n_a],"-m")==0)
            max_memory = atoi(argv[n_a+1]);

        if (strcmp(argv[n_a],"-d")==0)
            max_disk_space = atoi(argv[n_a+1]);

        if (strcmp(argv[n_a],"-v")==0)
            verbose = 1;

        if (strcmp(argv[n_a],"-vv")==0)
            verbose = 2;
        
        if (strcmp(argv[n_a],"-histo")==0)
            output_histo =true;
    }

    if (max_memory > total_ram)
    {
        printf("Maximum memory (%d MB), exceeds total RAM (%d MB). Setting maximum memory to %d MB.\n",max_memory,total_ram,total_ram/2);
        max_memory = total_ram/2;
    }

    STARTWALL(0);

    sorting_count(Reads,prefix,max_memory,max_disk_space,true,verbose);

    STOPWALL(0,"Total");

    delete Reads;

    return 0;
}


