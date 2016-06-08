output: main.o fastq_read.o
	g++ -std=c++11 main.o fastq_read.o -o umi
main.o: main.cpp
	g++ -c -std=c++11 main.cpp
fastq_read.o: fastq_read.cpp
	g++ -c -std=c++11 fastq_read.cpp
clean:
	rm *.o umi
