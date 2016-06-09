#ifndef _FASTQ_READ_H
#define _FASTQ_READ_H

class fastq_read {
public:
	/*
	  Default constructor
	*/
	fastq_read();

	/*
	  Creates a new fastq read

	  \param seq	Sequence read
	  \param plus	+
	  \param qual	Quality scores for read
	*/
	fastq_read(std::string seq, std::string plus, std::string qual);

	std::string _seq;
	std::string _plus;
	std::string _qual;
	int _count;
	bool _skip;
};

#endif // _FASTQ_READ_H
