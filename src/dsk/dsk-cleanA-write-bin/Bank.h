//
//  Bank.h
//
//  Created by Guillaume Rizk on 28/11/11.
//  Modified by Rayan Chikhi on 16/2/13
//

#ifndef Bank_h
#define Bank_h
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h> // Added by Pierre Peterlongo on 02/08/2012.

#define TAILLE_NOM 1024
#define MAX_NB_FILES 1000
#define BUFFER_SIZE 16384 // same as kseq.h 
#define WRITE_BUFFER  32768// 16384  //800000
#define KMERSBUFFER_MAX_READLEN 4096 // grows dynamically if needed
#define BINREADS_BUFFER 100000

using namespace std;

off_t fsize(const char *filename) ;

// heavily inspired by kseq.h from Heng Li (https://github.com/attractivechaos/klib)
typedef struct 
{
    gzFile stream;
    unsigned char *buffer;
    int buffer_start, buffer_end;
    bool eof;
    char last_char;
} buffered_file_t;

typedef struct 
{
    int length ,max;
    char *string;
} variable_string_t;
	
// supports opening multiple fasta/fastq files
class Bank{

    public:

        Bank(char *fname);
        Bank(char **fname, int nb_files_);
        void init(char **fname, int nb_files_);
        void close();

        bool get_next_seq(char **nseq, int *len);
        bool get_next_seq_from_file(char **nseq, int *len, int file_id);
    
        bool get_next_seq_from_file(char **nseq, char **cheader, int *len, int *hlen, int file_id);
        bool get_next_seq(char **nseq, char **cheader, int *len, int *hlen);

        void rewind_all();

        variable_string_t *read, *dummy, *header;

        int nb_files; // total nb of files
        int index_file; // index of current file
        uint64_t filesizes; // estimate of total size for all files

        signed int buffered_gets(buffered_file_t *bf, variable_string_t *s, char *dret, bool append, bool allow_spaces);

        ~Bank();

        uint64_t estimate_kmers_volume(int k);
        uint64_t estimate_nb_reads();
        int estimate_max_readlen();

        buffered_file_t  *buffered_file[MAX_NB_FILES];
};

class BinaryReads
{
    char filename[TAILLE_NOM];
    // const int sizeElement;
    unsigned char * buffer;
    int cpt_buffer;
    unsigned int  read_write_buffer_size;

    public:
    FILE * binary_read_file;

    BinaryReads(char *filename, bool write);
    // void write_element(void *element);
    //size_t read_element(void *element);
    void write_read(char * read, int readlen);
    void rewind_all();
    void close();
    void open(bool write);
    ~BinaryReads();

};

#endif
