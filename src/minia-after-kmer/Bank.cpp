//
//  Bank.cpp
//
//  Created by Guillaume Rizk on 28/11/11.
//

//TEST

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64

#include <algorithm>
#include <iostream>
#include <sys/stat.h>
#include <inttypes.h>
       #include <stdio.h>
       #include <stdlib.h>
       #include <string.h>
#include <cmath> // for log2f

#include "Bank.h"
#include "Kmer.h" // Bank (almost) doesn't need Kmer.h, but KmersBuffer certainly does
#include "lut.h"
#include <errno.h>
using namespace std;

off_t fsize(const char *filename) {
    struct stat st; 
    
    if (stat(filename, &st) == 0)
        return st.st_size;
    
    return -1; 
}

// the following functions are adapted from kseq.h by Heng Li (https://github.com/attractivechaos/klib)
inline bool rebuffer(buffered_file_t *bf)
{
    if (bf->eof)
        return false;
    bf->buffer_start = 0;
    bf->buffer_end = gzread(bf->stream, bf->buffer, BUFFER_SIZE);
    if (bf->buffer_end < BUFFER_SIZE)
        bf->eof = 1;
    if (bf->buffer_end == 0) 
        return false;
    return true;
}

inline signed char buffered_getc(buffered_file_t *bf)
{
    if (bf->buffer_start >= bf->buffer_end)
        if (! rebuffer(bf))
            return -1;
    return (signed char) ( bf->buffer[bf->buffer_start++] );
}

#define nearest_power_of_2(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

inline signed int Bank::buffered_gets(buffered_file_t *bf, variable_string_t *s, char *dret, bool append, bool allow_spaces)
{
    if (dret) *dret = 0;
    if (!append)
        s->length = 0;
    if (bf->buffer_start >= bf->buffer_end && bf->eof)
        return -1;
    while (1)
    {
        int i;
        if (bf->buffer_start >= bf->buffer_end)
            if (! rebuffer(bf))
                break;
        if (allow_spaces)
        {
            for (i = bf->buffer_start; i < bf->buffer_end ; i++)
                if (bf->buffer[i] == '\n') 
                    break;
        }
        else
        {
            for (i = bf->buffer_start; i < bf->buffer_end ; i++)
                // isspace() answers yes for ' ', \t, \n, \v, \f, \r
                if (isspace(bf->buffer[i]))
                    break;
        }
        if (s->max - s->length < (i - bf->buffer_start + 1))
        {
            s->max = s->length + (i - bf->buffer_start + 1);
            nearest_power_of_2(s->max);
            s->string = (char*)realloc(s->string,s->max);
        } 
        memcpy(s->string + s->length, bf->buffer + bf->buffer_start, i - bf->buffer_start);
        s->length += i - bf->buffer_start;
        bf->buffer_start = i + 1;
        if (i < bf->buffer_end)
        {
            if (dret)
                *dret = bf->buffer[i];
            break;
        }
    }
    if (s->string == NULL)
    {
        s->max = 256;
        s->string = (char*)calloc(256,1);
    }
    else if ( allow_spaces && s->length > 1 && s->string[s->length-1] == '\r')
        s->length--;
    s->string[s->length]= '\0';
    return s->length;
}

void Bank::rewind_all()
{
    for (int i=0; i<nb_files; i++)
    {
    	gzrewind(buffered_file[i]->stream);
        buffered_file[i]->last_char = buffered_file[i]->eof = buffered_file[i]->buffer_start = buffered_file[i]->buffer_end = 0;
    }
    index_file = 0;
}

// THIS READS FASTQ or FASTA, compressed with gzip or not
// no limit on read length, allows multi-line reads
// returns true if a read was successfuly read
//         false if end of file
// adapted from kseq.h by Heng Li (https://github.com/attractivechaos/klib)
bool  Bank::get_next_seq_from_file(char **nseq, char **cheader, int *len, int *hlen, int file_id)
{
    signed char c;
    buffered_file_t *bf = buffered_file[file_id];
    if (bf->last_char == 0)
    {
        while ( (c = buffered_getc(bf)) != -1 && c != '>' && c != '@'); // go to next header
        if (c == -1)
            return false; // eof
        bf->last_char = c;
    }
    read->length = dummy->length = 0;

    if (buffered_gets(bf, header, (char *)&c, false, false) < 0) //ici
        return false; // eof
    if (c != '\n')
        buffered_gets(bf, dummy, NULL, true, true); // read header //dummy instead of header to stop before first space
    
    if (read->string == NULL)
    {
        read->max = 256;
        read->string = (char*) malloc(read->max);
    }
    while ( (c = buffered_getc(bf)) != -1 && c != '>' && c != '+' && c != '@')
    {
        if (c == '\n')
            continue; // empty line
        read->string[read->length++] = c;
        buffered_gets(bf, read, NULL, true, true);
    }
    if (c == '>' || c == '@')
        bf->last_char = c;
    if (read->length + 1 >= read->max)
    {
        read->max = read->length + 2;
        nearest_power_of_2(read->max);
        read->string = (char*) realloc(read->string, read->max);
    }
    read->string[read->length] = '\0';
    if (c == '+') // fastq
    {
        if (dummy->max < read->max) // resize quality to match read length
        {
            dummy->max = read->max;
            dummy->string = (char*)realloc(dummy->string, dummy->max);
        }
        while ( (c = buffered_getc(bf)) != -1 && c != '\n'); // read rest of quality comment
        while (buffered_gets(bf, dummy, NULL, true, true) >= 0 && dummy->length < read->length); // read rest of quality
        bf->last_char = 0;
    }
    *len = read->length;
    *nseq = read->string;
    if (cheader && hlen)
    {
        *cheader = header->string;
        *hlen = header->length;
    }

    return true;
}

// wrapper
bool  Bank::get_next_seq_from_file(char **nseq, int *len, int file_id)
{
    return get_next_seq_from_file(nseq,NULL,len,NULL,file_id);
}

// wrapper
bool Bank::get_next_seq(char **nseq, char **cheader, int *len, int *hlen)
{
    bool success = get_next_seq_from_file(nseq,cheader,len,hlen,index_file);
    if (success)
        return true;
    
    // cycle to next file if possible
    if ( index_file < nb_files-1 )
    {
        index_file++;
        return get_next_seq(nseq,cheader, len,hlen);
    }
    return false;
}

// wrapper
bool Bank::get_next_seq(char **nseq, int *len)
{
  return get_next_seq(nseq,NULL,len,NULL);
}

// had to move the Bank(x,x) constructor to an init() to avoid calling a constructor inside the Bank(x) constructor
void Bank::init(char **fname, int nb_files_)
{
    int64_t i;
    nb_files = nb_files_;
    filesizes = 0;

    // open the reads file, don't know if it is a fasta/q file or a list of file names yet
    gzFile tempfile = gzopen(fname[0],"r");
    if (tempfile == NULL)
    {
        char *buffer = (char*)malloc(BUFSIZ);
        char * errorMessage = (char * ) strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s   %s \n",errorMessage,fname[0]);
        free(buffer);
        exit(1);
    }
    char deb=(char)gzgetc(tempfile);

    char **nfname;// [MAX_NB_FILES][TAILLE_NOM];
    nfname = (char**) malloc(sizeof(char*)*MAX_NB_FILES);
    for(int jj=0; jj<MAX_NB_FILES; jj++ )
	nfname [jj] =  (char*) malloc(sizeof(char)*TAILLE_NOM);
      
    if(deb=='>' || deb=='@' || deb==EOF)
    { // file is a fasta/q file
        gzclose(tempfile);
    }
    else // file contains a list of file names
    {
        char* ret;
        gzungetc(deb,tempfile);
        printf("File %s starts with character \"%c\", hence is interpreted as a list of file names\n",fname[0],deb );
        int ii;
        // get the filenames
        for (ii=0; ii<MAX_NB_FILES ; ii++)
        {
            ret = gzgets(tempfile, nfname[ii], BUFFER_SIZE);
            if (ret != NULL) {
                // remove \r \n chars
                char *endline = strchr(nfname[ii], '\n');
                if (endline)
                    *endline='\0';
                endline = strchr(nfname[ii], '\r');
                if (endline)
                    *endline='\0';
            }
            else // no more filenames
                break;
        }
        printf("Reading %i read files\n",ii);
        if(ii==MAX_NB_FILES)
            printf("Warning! using max number of read files (%i)\n",ii);

        nb_files = ii;
        fname = (char **) nfname;
        gzclose(tempfile);

    }

    // estimate total size of files
    for (i=0; i<nb_files; i++)
    {
        bool compressed = false;
        uint64_t estimated_filesize;

        if (strstr(fname[i],"gz") == (fname[i]+strlen(fname[i])-2) ) compressed=true;
        if (compressed)
            // crude hack, based on Quip paper reporting compression ratio (~0.3). 
            // gzseek(SEEK_END) isn't supported. need to read whole file otherwise :/
            estimated_filesize = fsize(fname[i]) * 4; 
        else
            estimated_filesize = fsize(fname[i]);

        filesizes += estimated_filesize;
    }

    // open each file for reading
    for (i=0; i<nb_files; i++)
    {
        buffered_file[i] = (buffered_file_t *)calloc(1, sizeof(buffered_file_t)); 
        buffered_file[i]->buffer = (unsigned char*) malloc(BUFFER_SIZE); 
        buffered_file[i]->stream = gzopen(fname[i],"r");

        if (buffered_file[i]->stream == NULL)
        {
            printf("error opening file: %s\n",fname[i]);
            exit(1);
        }
    }

    index_file = 0; // initialize the get_next_seq iterator to the first file

    // init read and dummy (for readname and quality)
    read = (variable_string_t*) calloc(1,sizeof(variable_string_t));
    dummy = (variable_string_t*) calloc(1,sizeof(variable_string_t));
    header = (variable_string_t*) calloc(1,sizeof(variable_string_t));

    
    for(int jj=0; jj<MAX_NB_FILES; jj++ )
      free	(nfname [jj]); 
    free(nfname);
}

Bank::Bank(char *fname0)
{
    char *fname[1] = { fname0 };
    init(fname, 1);
}


Bank::Bank(char **fname, int nb_files_)
{
    init(fname,nb_files_);
}

Bank::~Bank(){
    variable_string_t * to_free[3] = {read, dummy, header};
    for (int i = 0; i < 3; i++)
    {
        if (to_free[i])

        {
            if (to_free[i]->string)
                free(to_free[i]->string);
            free(to_free[i]);
        }
    }
    for (int i=0; i<nb_files; i++)
    {
        free(buffered_file[i]->buffer);
        free(buffered_file[i]);
    }
}

void Bank::close()
{
    for (int i=0; i<nb_files; i++)
        gzclose(buffered_file[i]->stream);
}

// estimate the volume of all redundant kmers in the reads, if they were to be stored in 2bits
uint64_t Bank::estimate_kmers_volume(int k)
{

    char * rseq;
    int readlen;
    int NbRead = 0;
    //int kmer_nbits = std::max(64,(int)pow(2,ceilf(log2f(2*k)))); // Bank assumes that a kmer is stored in the smallest integer type (e.g. uint64_t or uint128_t) // not accurate anymore with _ttmath/_largeint
    int kmer_nbits = sizeof(kmer_type)*8;
    rewind_all();
    uint64_t volume = 0;

    while (get_next_seq(&rseq,&readlen))
    {
        if (readlen >= k)
            volume += (readlen-k+1) * (uint64_t) kmer_nbits;
        if (NbRead++ == 1000)
            break;
    }

    if ( gztell(buffered_file[index_file]->stream) == 0) // empty file
        return 1;

    volume = volume * ((float)filesizes/gztell(buffered_file[index_file]->stream));

    volume = volume / 1024 /1024 /8; // put it in MB
    
    if (volume == 0)  // tiny files fix
        volume = 1;

    rewind_all();
    return volume;
}

// estimate the number of reads
uint64_t Bank::estimate_nb_reads()
{
    char * rseq;
    int readlen;
    int NbRead = 0;
    rewind_all();
    
    uint64_t volume = 0;
    while (get_next_seq(&rseq,&readlen))
    {
        volume += 1;
        if (NbRead++ == 1000)
            break;
    }

    if ( gztell(buffered_file[index_file]->stream) == 0) // empty file
        return 1;

    volume = (volume * filesizes) / gztell(buffered_file[index_file]->stream); // linear extrapolation from the first 1k reads 

    rewind_all();
    return volume;
}

// estimate maximum read length
// from the first 10000 reads of each file
int Bank::estimate_max_readlen()
{
    char * rseq;
    int readlen;
    rewind_all();
    int max_readlen = 0;
    uint64_t volume = 0;

    index_file = 0;
    
    while ( index_file < nb_files )
    {
        int NbRead = 0;
        while (get_next_seq_from_file(&rseq,NULL,&readlen,NULL,index_file))
        {
            max_readlen = max(readlen, max_readlen);
            if (NbRead++ == 10000)
                break;
        }
        index_file++;
    } 

    index_file = 0;
    return max_readlen;
}

// BinaryBank: a binary file containing kmers

BinaryBank::BinaryBank(char *given_filename, int given_sizeElement, bool write) : sizeElement(given_sizeElement)
{
    strcpy(filename,given_filename);
    open(write);
    buffer_size_nelem= (WRITE_BUFFER/given_sizeElement);
    buffer = (void *) malloc(given_sizeElement * buffer_size_nelem);
    cpt_buffer=0;
}


BinaryBankConcurrent::BinaryBankConcurrent(char *given_filename, int given_sizeElement, bool write, int given_nthreads) : BinaryBank(given_filename,given_sizeElement,write) 
{
    nthreads = given_nthreads;
    
    //free(buffer); buffer =NULL; //cannot do that
    bufferT = (void **) malloc(sizeof(void*) * nthreads);

    for (int i= 0; i< nthreads; i++)
    {
         ((void ** )bufferT)[i]= (void *) malloc( WRITE_BUFFER);
      //  ((void ** )bufferT)[i]= (void *) malloc(sizeElement* WRITE_BUFFER);

    }
    cpt_buffer_tid = (int  *)malloc(sizeof(int) * nthreads);
    memset (cpt_buffer_tid,0,sizeof(int) * nthreads);
}


void BinaryBankConcurrent::write_element_buffered( void *element, int tid)
{
    
    if(cpt_buffer_tid[tid]>= WRITE_BUFFER -100)
    {
        flush(tid);
    }
    
   // ((kmer_type **)bufferT)[tid][ cpt_buffer_tid[tid] / sizeElement]= *((kmer_type *)element);
    
    
    char * buf_pt = ((char**) bufferT)[tid];
    memcpy(buf_pt + cpt_buffer_tid[tid] , element, sizeElement);
    cpt_buffer_tid[tid]+=sizeElement;
    
    
//    char * buf_pt = ((char**) bufferT)[tid];
//    buf_pt +=  cpt_buffer_tid[tid];
//    kmer_type * write_adress =  (kmer_type *) buf_pt;
//    *write_adress = *((kmer_type *)element);
//    
//   //   *((kmer_type *) (&(((char **)bufferT) [tid] [cpt_buffer_tid[tid]] ))) = *((kmer_type *)element); //works but ugly and may break
//
//    cpt_buffer_tid[tid]+=sizeElement;

    
}


void BinaryBankConcurrent::write_buffered( void *element, int size, int tid)
{
    write_buffered( element, size, tid, true);
}

void BinaryBankConcurrent::write_buffered( void *element, int size, int tid, bool can_flush)
{
    if(cpt_buffer_tid[tid]>= WRITE_BUFFER -100 && can_flush)
    {
        flush(tid);
    }
    
    char * buf_pt = ((char**) bufferT)[tid];    
    memcpy(buf_pt + cpt_buffer_tid[tid] , element, size);
    
    cpt_buffer_tid[tid]+=size;
    // cpt_buffer_tid[tid]++;
    

}



void BinaryBankConcurrent::flush(int tid)
{
    flockfile(binary_read_file);
    if (!fwrite( ((void **)bufferT)[tid], 1, cpt_buffer_tid[tid], binary_read_file))            
    {
        printf("error: can't fwrite (disk full?)\n");
        funlockfile(binary_read_file);
        exit(1);
    }
    cpt_buffer_tid[tid]=0;
    funlockfile(binary_read_file);

}


//should be called by only one of the threads
void BinaryBankConcurrent::close()
{
    //flush buffer // if close Bank in read mode with data in the readbuffer, will result in error
    for(int ii=0; ii< nthreads; ii++)
    {
        if(cpt_buffer_tid[ii])
        {
            if (!fwrite(((void **)bufferT)[ii], 1, cpt_buffer_tid[ii], binary_read_file))
          //      if (!fwrite(((void **)bufferT)[ii], sizeElement, cpt_buffer_tid[ii], binary_read_file))

            {
                printf("error: can't fwrite (disk full?)\n");
                exit(1);
            }
        }
        cpt_buffer_tid[ii]=0;
    }
    
    fclose(binary_read_file);
}


void BinaryBank::write_element( void *element)
{
  //  flockfile(binary_read_file);
   // fprintf(stderr,"write elem %lli \n",*(int64_t *)element);
    if (!fwrite(element, sizeElement, 1, binary_read_file))
    {
       // funlockfile(binary_read_file);
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
  //  funlockfile(binary_read_file);
}


void BinaryBank::write_element_buffered( void *element)
{
    
    if(cpt_buffer==buffer_size_nelem)
    {
        if (!fwrite(buffer, sizeElement, buffer_size_nelem, binary_read_file))
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
        cpt_buffer=0;
    }
    
    
    ((kmer_type *)buffer)[cpt_buffer]= *((kmer_type *)element);
    cpt_buffer++;
    
}



size_t BinaryBank::read_element( void *element)
{
    return fread(element, sizeElement,1, binary_read_file);
}

size_t BinaryBank::read_element_buffered( void *element)
{
    if(cpt_buffer==0)
    {
        cpt_buffer=fread(buffer, sizeElement,buffer_size_nelem, binary_read_file);
        if (cpt_buffer==0) return 0;
    }
    *((kmer_type *)element) =  ((kmer_type *)buffer)[cpt_buffer-1] ; // todo check read order is consisten with file
    cpt_buffer --;
    return cpt_buffer+1; // nb remaining before read
}

// used to read/write raw information to the binary file (e.g. kmer count)

void BinaryBank::write( void *element, int size)
{
    if (!fwrite(element, size, 1, binary_read_file))
    {
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
}

size_t BinaryBank::read( void *element, int size)
{
    return fread(element, size,1, binary_read_file);
}


void BinaryBank::rewind_all()
{
    rewind(binary_read_file);
}

void BinaryBank::close()
{
    //flush buffer // if close Bank in read mode with data in the readbuffer, will result in error
    if(cpt_buffer)
    {
    if (!fwrite(buffer, sizeElement, cpt_buffer, binary_read_file))
    {
        printf("error: can't fwrite (disk full?)\n");
        exit(1);
    }
    }
    cpt_buffer=0;
    
    fclose(binary_read_file);
}

void BinaryBank::open(bool write)
{
    binary_read_file = fopen(filename,write?"wb":"rb");
    if( binary_read_file == NULL )
    {
        char *buffer = (char*)malloc(BUFSIZ);
        char * errorMessage = (char * ) strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s  write %i  %s\n",errorMessage,write,filename);
        free(buffer);
        exit(1);
    }

}

off_t BinaryBank::nb_elements()
{
  return fsize(filename)/sizeElement;
}


BinaryBank::~BinaryBank()
{
    if(buffer!=NULL)
    {
        free (buffer); //buffer =NULL;
    }
}


BinaryBankConcurrent::~BinaryBankConcurrent()
{
    
    for (int i= 0; i< nthreads; i++)
    {
        free(((void ** )bufferT)[i]);
        ((void ** )bufferT)[i]=NULL;
    }
    free(bufferT);
}



/////////////class BinaryReads a file containing reads

BinaryReads::~BinaryReads()
{
    free (buffer); buffer = NULL;
}


BinaryReads::BinaryReads(char *given_filename,  bool write)
{
    read_write_buffer_size = BINREADS_BUFFER;
    strcpy(filename,given_filename);
    open(write);
    buffer = (unsigned char *) malloc(read_write_buffer_size*sizeof(unsigned char));
    cpt_buffer = 0;
}


void BinaryReads::rewind_all()
{
    rewind(binary_read_file);
}

void BinaryReads::close()
{
    unsigned int block_size =0;
    //flush buffer
    if(cpt_buffer)
    {
        //printf("close :write block %i \n",cpt_buffer);
        block_size = cpt_buffer;
        fwrite(&block_size, sizeof(unsigned int), 1, binary_read_file); // block header
        if (!fwrite(buffer, 1, cpt_buffer, binary_read_file))
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
    }
    cpt_buffer=0;
    
    fclose(binary_read_file);
}

void BinaryReads::open(bool write)
{
    binary_read_file = fopen(filename,write?"wb":"rb");
    if( binary_read_file == NULL )
    {
        char *buffer = (char*)malloc(BUFSIZ);
        char * errorMessage = (char * ) strerror_r( errno, buffer, BUFSIZ ); // get string message from errno
        printf("error during fopen: %s  write %i  %s\n",errorMessage,write,filename);
        free(buffer);
        exit(1);
    }
    
}



//format is
// 32 bit integer = readlen,  then seq in binary
// then next read..
//32 bit len is overkill but simpler
//also makes buffer then write block with header : size of block to read, with n reads .... will allow large fread when reading this file ...
void BinaryReads::write_read(char * read, int readlen)
{
    int tai = readlen;
    unsigned char rbin;
    char * pt = read;
    unsigned int block_size = 0;
    
 //   printf("write read %i / %i   readlen %i \n",cpt_buffer,read_write_buffer_size,readlen);
    //todo : also flush to disk  sometimes (ie if very large buffer, to create smaller blocks..)
    if(cpt_buffer >= (read_write_buffer_size-readlen) || cpt_buffer > 10000000 )  ////not enough space to store next read   true space is 4 + readlen/4 + rem
        //flush buffer to disk
    {
        
        block_size = cpt_buffer;
        
        //printf("write block %i\n",block_size);
        if(block_size) fwrite(&block_size, sizeof(unsigned int), 1, binary_read_file); // block header
        if (!fwrite(buffer, 1, cpt_buffer, binary_read_file)) // write a block, it ends at end of a read
        {
            printf("error: can't fwrite (disk full?)\n");
            exit(1);
        }
        cpt_buffer=0;
    }
    
    //check if still not enough space in empty buffer : can happen if large read, then enlarge buffer
    if(read_write_buffer_size < readlen)
    {
        read_write_buffer_size = 2*readlen; // too large but ok
        buffer =  (unsigned char *) realloc(buffer,sizeof(unsigned char) * read_write_buffer_size);
    }
    
    memcpy(buffer+cpt_buffer,&readlen,sizeof(int));
    cpt_buffer+= sizeof(int);
    
    //fwrite( (void *) &readlen, sizeof(int), 1, binary_read_file);

    
    for (tai=readlen; tai>=4  ; tai-=4)
    {
        rbin = code4NT(pt);
      //  fwrite((void *) &rbin, 1,1,binary_read_file );
        buffer[cpt_buffer]=rbin; cpt_buffer++;
        pt +=4;
    }
    
    //then remaining
    if(tai)
    {
        rbin = code_n_NT(pt,tai);
       // fwrite( (void *) &rbin,1,1,binary_read_file);
        buffer[cpt_buffer]=rbin; cpt_buffer++;
    }
}



void  compute_kmer_table_from_one_seq(int readlen, char * seq, kmer_type * kmer_table )  //,char * pkmer_table //pour remplissage table loc
{
    kmer_type graine = codeSeed(seq);
    kmer_type graine_revcomp = revcomp(graine);
    kmer_table[0] = min(graine,graine_revcomp);
    seq++;
    for (int i=1; i<readlen-sizeKmer+1; i++)
    {
        graine =   (graine * 4 + NT2int(seq[sizeKmer-1])) & kmerMask   ;
        graine_revcomp =  ((graine_revcomp >> 2) +  ( ((kmer_type) comp_NT[NT2int(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask ;
        kmer_table[i] = min(graine,graine_revcomp);
        seq++;
    }
}





////kmers buffer


KmersBuffer::KmersBuffer(BinaryReads *bfile, int  pbuffer_size, int nseq_task )
{
    
    read_write_buffer_size = BINREADS_BUFFER;
    buffer = ( char *) malloc(read_write_buffer_size*sizeof( char));
    cpt_buffer = 0;
    cpt_binSeq_read =0; binSeq_toread =0;
    max_read_length = KMERSBUFFER_MAX_READLEN;
    binfile = bfile;
    buffer_size = pbuffer_size;
    kmers_buffer =  (kmer_type *) malloc(sizeof(kmer_type) * buffer_size);
   // binSeq =  (char *) malloc(sizeof(char) * max_read_length); // no need to alloc ram for binse : will points to buffer
    binSeq_extended =  (char *) malloc(sizeof(char) * max_read_length);
    blocksize_toread =0;


    nseq_step = nseq_task;
    binary_read_file = bfile->binary_read_file;
    
}


void KmersBuffer::reset_max_readlen(int read_length)
{

    max_read_length = read_length;

  //  binSeq =  (char *) realloc(binSeq,sizeof(char) * max_read_length);
    binSeq_extended =  (char *) realloc(binSeq_extended,sizeof(char) * max_read_length);

}

 KmersBuffer::~KmersBuffer()
{
    free (kmers_buffer);
    free(buffer);
    //free(binSeq);
    free(binSeq_extended);
}


//now returns number of kmers read
int KmersBuffer::readkmers()
{
    
//    printf("--------------\ncalling readkmers\n");
//    printf("cpt_buffer %i blocksize_toread %i\n",cpt_buffer,blocksize_toread);
//    printf("cpt_binSeq_read %i binSeq_toread %i\n------------\n",cpt_binSeq_read,binSeq_toread);

    int llen;
    int * len = & llen ;
    unsigned int block_size =0;
    
    //////reading new block from disk if needed
    if(cpt_buffer == blocksize_toread  && (binSeq_toread  <= cpt_binSeq_read))
    {
        flockfile(binary_read_file);
        if( ! fread(&block_size,sizeof(unsigned int),1, binary_read_file)) //read block header
        {
            funlockfile(binary_read_file);
       //     printf("trying to read new block but end : return 0\n");
            return 0; // no more blocks to read
        }

        
        if(block_size >= read_write_buffer_size)
        {
            read_write_buffer_size = 2*block_size;
            buffer =  ( char *) realloc(buffer,sizeof( char) * read_write_buffer_size);
        }
        
        fread(buffer,sizeof( char),block_size, binary_read_file); // read a block of reads into the buffer
        funlockfile(binary_read_file);
        cpt_buffer = 0;
        blocksize_toread = block_size;
      //  printf("reading block %i  %i/%i\n",block_size,cpt_buffer, blocksize_toread);

    }
    ///////////////////////
    
    //now parse the whole block in ram
    
    int i,j;
    int nchar;
    unsigned char fournt;
    
    nkmers = 0;
    int nseq_lues = 0;
    while(cpt_buffer < blocksize_toread || ( binSeq_toread  > cpt_binSeq_read)) //while work to do
    {
       // printf("cpt_buffer %i blocksize_toread %i\n",cpt_buffer,blocksize_toread);

        if( binSeq_toread <= cpt_binSeq_read)// read new read if needed
        {
            memcpy(len,buffer+cpt_buffer,sizeof(int)); // read len
            cpt_buffer += sizeof(int);
            nseq_lues ++;
            
            if( (*len) > max_read_length) reset_max_readlen((int)(1.2*(*len))); // in ram 2 times the max size of reads
            nchar = ((*len)+3)/4;
            // fread(binSeq, sizeof(char), nchar, binary_read_file ); // read one seq from binfile
            binSeq = buffer + cpt_buffer; // point binseq to correct place
            cpt_buffer += nchar;
            
            j=0;
            for(i=0; i<nchar; i++)
            {
                fournt = binSeq[i];
                binSeq_extended[j+3]=fournt & 3; fournt = fournt >> 2; // il faudrait deporter ce calcul du lock .. cest fait
                binSeq_extended[j+2]=fournt & 3; fournt = fournt >> 2;
                binSeq_extended[j+1]=fournt & 3; fournt = fournt >> 2;
                binSeq_extended[j+0]=fournt & 3;
                j+=4;
            }
            binSeq_toread = *len-sizeKmer+1;
            cpt_binSeq_read = 0;
            
//            printf("Init cpt_binSeq_read %i  binSeq_toread %i  todo ? %i\n",cpt_binSeq_read,binSeq_toread,*len-sizeKmer+1);
//            printf("cpt_buffer %i blocksize_toread %i\n",cpt_buffer,blocksize_toread);

        }
        
        
        {
   //         printf("cpt_binSeq_read %i binSeq_toread %i\n",cpt_binSeq_read,binSeq_toread);

            //printf("cpt_binSeq_read %i  binSeq_toread %i \n",cpt_binSeq_read,binSeq_toread);
            char *seq = binSeq_extended+cpt_binSeq_read;
            kmer_type graine;
            kmer_type graine_revcomp;
            if( binSeq_toread  > cpt_binSeq_read)
            {
                graine = codeSeed_bin(seq);
                graine_revcomp = revcomp(graine);
            
            if(nkmers>=buffer_size)
            {
                return nkmers;
            }
            kmers_buffer[nkmers] = min(graine,graine_revcomp); nkmers++; cpt_binSeq_read ++;
            seq++;
            }
            
            while( binSeq_toread  > cpt_binSeq_read)
            {
               // printf(" here cpt_binSeq_read %i binSeq_toread %i\n",cpt_binSeq_read,binSeq_toread);

                graine =  (graine * 4 + (seq[sizeKmer-1])) & kmerMask ;
                graine_revcomp =  ((graine_revcomp >> 2) +  ( ((kmer_type) comp_NT[(int)(seq[sizeKmer-1])]) <<  (2*(sizeKmer-1))  )  ) & kmerMask;
                kmers_buffer[nkmers] = min(graine,graine_revcomp); nkmers ++; cpt_binSeq_read ++;
                seq++;
                if(nkmers>=buffer_size)
                {
                  //  printf(" return ... %i / %i \n",nkmers,buffer_size);

                    return nkmers;
                }
            }
            
            //return 1; //  en fait lecture partielle d'un read
            
        }
        
        
    }
    
  //  printf("return nseq lues %i  nkmers %i\n",nseq_lues,nkmers);

    return nkmers;
    
    
}

 
