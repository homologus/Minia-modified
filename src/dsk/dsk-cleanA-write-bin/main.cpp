#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include "Bank.h"

using namespace std;

char filename[TAILLE_NOM];
unsigned char * buffer;
int cpt_buffer;
int max_read_length = KMERSBUFFER_MAX_READLEN;
int readlen;
const char *binary_read_file = (char *)"reads_binary";

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


int main(int argc, char *argv[])
{
	Bank *Sequences = new Bank(argv[1]);

    // default prefix is the reads file basename
    char *reads_path=strdup(argv[1]);
    string reads_name(basename(reads_path)); // posix basename() may alter reads_path
    free(reads_path);
    int lastindex = reads_name.find_last_of(".");
    strcpy(prefix,reads_name.substr(0, lastindex).c_str());


	BinaryReads *binread = NULL;
	binread = new BinaryReads(return_file_name(binary_read_file),true);

	//start by the conversion of the file to binary format

	char * pt_begin;
	int idx =0 ;
	int64_t NbRead = 0;
	char * rseq;

	Sequences->rewind_all();
	while(1)
	{
		if(! Sequences->get_next_seq(&rseq,&readlen)) break; // read  original fasta file

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

		NbRead++;
	}
	binread->close();

	return 0;
} 

