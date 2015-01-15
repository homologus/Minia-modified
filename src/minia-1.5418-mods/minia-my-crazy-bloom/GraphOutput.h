#include <tr1/unordered_map>
#include <tr1/functional>

#include <set>
#include <stdlib.h> // for exit()
#include <iostream>
#include <string>
#include <fstream>
#include <regex.h>
#include <stdio.h>
#include <assert.h> 

#include "Kmer.h"
#include "Bank.h"

#ifndef _GRAPHOUTPUT_H
#define _GRAPHOUTPUT_H

using namespace std;
using namespace tr1;

// hash functions for unordered_map with various kmer_type's
namespace std
{

   //structure for print id nodes and edges in graph output
   struct id_els{
		long node;
		long edge;
	       };
namespace tr1
{

    #ifdef _LP64
   template <>
   struct hash<__uint128_t> : public unary_function<__uint128_t, size_t>
   {
       size_t operator()(const __uint128_t& elem) const
       {
           hash<uint64_t> hash_func;
           return hash_func((uint64_t)(elem>>64)) ^ hash_func((uint64_t)(elem&((((__uint128_t)1)<<64)-1)));
       }
   };
    #endif

    #ifdef _ttmath
   template <>
   struct hash<ttmath::UInt<KMER_PRECISION> > : public unary_function<ttmath::UInt<KMER_PRECISION>, size_t>
   {
       size_t operator()(const ttmath::UInt<KMER_PRECISION>& elem) const
       {
           hash<uint64_t> hash_func;

           // hash = XOR_of_series[hash(i-th chunk iof 64 bits)
           uint64_t result = 0, to_hash;
           ttmath::UInt<KMER_PRECISION> intermediate = elem;
           uint32_t mask=~0, chunk;
           int i;
           for (i=0;i<KMER_PRECISION/2;i++)
           {
               // retrieve a 64 bits part to hash 
               (intermediate & mask).ToInt(chunk);
               to_hash = chunk;
               intermediate >>= 32;
               (intermediate & mask).ToInt(chunk);
               to_hash |= ((uint64_t)chunk) << 32 ;
               intermediate >>= 32;

               result ^= hash_func(to_hash);
           }
           return result;
       }
   };
    #endif

    #ifdef _largeint
   template <>
       struct hash<LargeInt<KMER_PRECISION> > : public unary_function<LargeInt<KMER_PRECISION>, size_t>
       {
           size_t operator()(const LargeInt<KMER_PRECISION>& elem) const
           {
               hash<uint64_t> hash_func;

               // hash = XOR_of_series[hash(i-th chunk iof 64 bits)
               uint64_t result = 0, to_hash;
               LargeInt<KMER_PRECISION> intermediate = elem;
               uint32_t mask=~0, chunk;
               int i;
               for (i=0;i<KMER_PRECISION/2;i++)
               {
                   // retrieve a 64 bits part to hash 
                   chunk = (intermediate & mask).toInt();
                   to_hash = chunk;
                   intermediate = intermediate >> 32;
                   chunk = (intermediate & mask).toInt();
                   to_hash |= ((uint64_t)chunk) << 32 ;
                   intermediate = intermediate >> 32;

                   result ^= hash_func(to_hash,num_hash);
               }
               return result;
           }
       };
    #endif
}
}

class GraphOutput {

public:
    string prefix;
    string graph_file_name;
    string nodes_file_name;
    string edges_file_name;
    string xml_file_name;
    string json_nodes_file_name;
    string json_edges_file_name;
    string json_file_name;
    int graph_format;
    id_els first_id_els, el_graph;

    long edge_id; // the json format needs an id on the nodes. 
  
    static const string graph_file_suffix;
    static const string nodes_file_suffix;
    static const string edges_file_suffix;
    static const string xml_file_suffix;
    static const string json_nodes_file_suffix;
    static const string json_edges_file_suffix;
    static const string json_file_suffix;
    
    FILE *graph_file,*nodes_file,*edges_file;


    GraphOutput(string prefix, int graph_format);
    GraphOutput(string prefix, int graph_format, id_els first_id_els); //PIERRE
    void close();

    long sequence_length(string line);
    int print_node(long index, char *ascii_nodet);
    void print_edge(long index, long id, long id2, string label);


    enum LeftOrRight { LEFT=0, RIGHT=1 };
    enum Strand { FW=0, RC=1 };
    struct node_strand {
        long node;
        Strand strand;
        LeftOrRight left_or_right;
        node_strand(long node, Strand strand, LeftOrRight left_or_right) : node(node), strand(strand), left_or_right(left_or_right) {}
        bool operator<(const node_strand &other) const {
            if (node != other.node)
                return (node < other.node);
            if (left_or_right != other.left_or_right)
                return left_or_right < other.left_or_right;
            return (strand < other.strand);
        }
    };

    std::tr1::unordered_map<kmer_type,set<node_strand> > kmer_links;

    id_els construct_graph(string linear_seqs_name); // PIERRE: added the return value
    void load_nodes_extremities(string linear_seqs_name);

 private: 
    void init(bool erase); // PIERRE
};
#endif //_GRAPHOUTPUT_H
