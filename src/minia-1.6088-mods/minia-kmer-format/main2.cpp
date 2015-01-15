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

#include "Bank.h"

int main()
{
	kmer_type kmer;
	char seq[46];
	char line[128];
	size_t len = 0;

	FILE *fpt=fopen("kmers","r");
	BinaryBankConcurrent * SolidKmers = new BinaryBankConcurrent("save",sizeof(kmer),true,1);

	while ( fgets ( line, 1000, fpt ) != NULL ) 
	{
		sscanf(line,"%s",seq);
		kmer=codeSeed(seq);
		SolidKmers->write_element_buffered(&kmer,0);
	}
	SolidKmers->close();
}
