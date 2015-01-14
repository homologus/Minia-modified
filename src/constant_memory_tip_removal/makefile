CC=g++
PROGNAME=constant_memory_tip_removal
CFLAGS =  -O4 -DNO_BLOOM_UTILS
SRC=minia/Bank.cpp minia/Kmer.cpp minia/OAHash.cpp minia/Utils.cpp
EXEC=$(PROGNAME)
OBJ= $(SRC:.cpp=.o)

ifeq ($(prof),1)
 CFLAGS+=-pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g 
endif

ifdef k
 KMER_SIZE := $(k)
else
 KMER_SIZE := 64
endif


K_BELOW_32 := $(shell echo $(KMER_SIZE)\<=32 | bc)
K_BELOW_64 := $(shell echo $(KMER_SIZE)\<=64 | bc)
ARCH := $(shell getconf LONG_BIT) # detects sizeof(int)
USING_UINT128 := 0
ttmath := 0

ifeq ($(K_BELOW_32),0)

    # use uint128 when k<=64 and 64-bit architecture
    ifeq ($(K_BELOW_64),1)
        ifeq ($(strip $(ARCH)),64)
            CFLAGS += -Dkmer_type=__uint128_t
            USING_UINT128 := 1
        endif
    endif
    
    # use a bigint library otherwise
    ifeq ($(USING_UINT128),0)
        ttmath := 1
    endif
endif

# ttmath is used when you type "make k=[kmer size]" with a kmer size longer than supported integer type,
# or when typing "make k=[anything] ttmath=1"
ifeq ($(ttmath),1)
    KMER_PRECISION := $(shell echo \($(KMER_SIZE)+15\)/16 | bc)
    CFLAGS += -D_ttmath -DKMER_PRECISION=$(KMER_PRECISION)
endif

all: $(EXEC) 

$(PROGNAME): clean $(OBJ) main.cpp
	$(CC) -o $@ $(OBJ) main.cpp $(CFLAGS) -lz

%.o: %.cpp %.h
	$(CC) -o $@ -c $< $(CFLAGS)


%.o: %.c %.h 
	$(CC) -o $@ -c $< $(CFLAGS)
    
clean:
	rm -rf *.o minia/*.o

test: all
	cd simple_test && ../constant_memory_tip_removal lambdaphage_nodes.fasta lambdaphage.solid_kmers_binary_with_count 10
