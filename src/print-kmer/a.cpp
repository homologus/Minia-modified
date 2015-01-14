#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <inttypes.h>
#include <cmath> // for log2f
#include <algorithm> // for max
#include <unistd.h> // for truncate
#include "Kmer.h"

typedef uint64_t kmer_type;

main()
{
	kmer_type kmer;
	int i;
	char *S;
	sizeKmer=21;
	kmerMask=(((kmer_type)1)<<(sizeKmer*2))-1;
	FILE *T2_file = fopen("a.solid_kmers_binary", "r");
	fread(&i, 4, 1, T2_file);
	printf("%i\n",i);


	while (fread(&kmer, sizeof(kmer), 1, T2_file))
	{
		S=print_kmer(kmer,sizeKmer,kmerMask);
		printf("%s\n",S);
	}
}
