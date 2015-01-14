//
//  Bank.cpp
//
//  Created by Guillaume Rizk on 28/11/11.
//


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
#include <errno.h>
using namespace std;

int NT2int(char nt)
{
    int i;
    i = nt;
    i = (i>>1)&3; // that's quite clever, guillaume.
    return i;
}

unsigned char  code4NT(char *seq)
{
    int i;
    unsigned char x;

    x=0;
    for (i=0; i<4; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    return x;
}


unsigned char  code_n_NT(char *seq, int nb)
{
    int i;
    unsigned char x;

    x=0;
    for (i=0; i<nb; ++i)
    {
        x = x*4 + NT2int(seq[i]);
    }
    x = x << ((4-nb)*2) ;
    return x;
}


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


