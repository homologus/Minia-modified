/* The MIT License

   Copyright (c) 2013 Mohua Bose

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

*/


/* Contact: Mohua Bose <bose@homolog.us> */


#include <stdio.h>
#include <stdlib.h>
#include "Bank.h"

float NBITS_PER_KMER = 11 ; // number of bits per kmer that optimizes bloo1 size



//
//	main program
//

int main(int argc, char **argv)
{
	BinaryBank * SolidKmers;
	printf("%i\n",sizeof(kmer_type));
	SolidKmers = new BinaryBank("a.solid_kmers_binary",sizeof(kmer_type),0);
	printf("nelem %lli nbits %g \n",SolidKmers->nb_elements(),NBITS_PER_KMER);
}
