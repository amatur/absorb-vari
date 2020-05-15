CC=g++
CFLAGS=-c -std=c++11 -O2

SOURCES=decoder.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=decoder


BCALMFILE=/Volumes/exFAT/data2019/staphsub/31/list_reads.unitigs.fa
K=31
BCALMFILE=/Users/Sherlock/amaturWS/data2019/staphsub/$K/list_reads.unitigs.fa


all: $(SOURCES) $(EXECUTABLE)
	rm decoder.o

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

#rule to build sources
.cpp.o: *.hpp
	$(CC) $(CFLAGS) $< -o $@

#tell that any header change should cause a rebuild
#decoder.cpp: decoder.hpp


decoderd: decoder.cpp
	$(CC) $(CFLAGS) decoder.cpp -o decoder.o
	$(CC) decoder.o -o decoder.out
	
clearop:
	rm -f absorbDecompressed.fa tip.fa

clean:
	rm -f *.o decoder

clear:
	rm -f absorbDecompressed.fa tip.fa
	rm -f *.o decoder

.SILENT:run

test:
	./decoder -i $(PWD) -k 11 -t 0

