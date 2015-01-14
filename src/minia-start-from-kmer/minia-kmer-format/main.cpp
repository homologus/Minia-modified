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
	BinaryBank *SolidKmers = new BinaryBank("o.solid_kmers_binary",sizeof(kmer_type),0);
	//BinaryBank *SolidKmers = new BinaryBank("save",sizeof(kmer_type),0);
        char seq[22];

	SolidKmers->rewind_all();
	while (SolidKmers->read_element(&kmer))
	{
		code2seq (kmer, seq);
		printf("%s\n",seq);
	}
}

