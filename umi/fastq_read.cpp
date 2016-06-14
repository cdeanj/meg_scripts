#include <string>
#include "fastq_read.h"

fastq_read::fastq_read() {}

fastq_read::fastq_read(const std::string seq, const std::string plus, std::string qual) {
	_seq = seq;
	_plus = plus;
	_qual = qual;
	_count = 1;
}
