#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <libgen.h>
#include <string.h>
#include <map>
#include "args.h"
#include "fastq_read.h"
using namespace std;

string base_name(string fp) {
	return string(basename(&fp[0]));
}

/*
 * Passes a map of fastq reads and postfix datum. If the skip boolean for each
 * umi is set, then there must be a collection of fastq reads with identical umis
 * that differed in their sequences. These are written to a mismatch file.
 * UMI's with identical sequences are written to a match file
 */
void write_to_file(const map<string,fastq_read> &mapper, const string &postfix) {
	ofstream ofs_m("match_" + postfix);
	ofstream ofs_mm("mismatch_" + postfix);
	map<string,fastq_read>::const_iterator it;
	for(it = mapper.begin(); it != mapper.end(); ++it) {
		if(!it->second._skip) {
			ofs_m << "@" << it->first << "/" << it->second._count << endl;
			ofs_m << it->second._seq << endl;
			ofs_m << it->second._plus << endl;
			ofs_m << it->second._qual << endl;
		}
		else {
			ofs_mm << "@" << it->first << endl;
			ofs_mm << it->second._seq << endl;
			ofs_mm << it->second._plus << endl;
			ofs_mm << it->second._qual << endl;
		}
	}
	ofs_m.close();
	ofs_mm.close();
}

/*
 * Returns a 24-mer extracted from each each sequence identifier
*/
string get_umi(const string &line) {
	int idx = line.find("|")+1;
	return line.substr(idx,24);
}

/*
 * Returns a map of fastq reads, with a umi representing each sequence
 */
map<string, fastq_read> process_fastq(const string &f) {
	map<string, fastq_read> m;
	string line, umi, seq, plus, qual;

	ifstream in(f);
	while(getline(in, line)) {
		if(line[0] == '@') {
			umi = get_umi(line);
			getline(in, seq);
			getline(in, plus);
			getline(in, qual);	
		}
		if(m.count(umi) > 0 && (m[umi]._skip == false) ) {
			if(m[umi]._seq == seq) {
				m[umi]._count += 1;
			}
			else {
				m[umi]._skip = true;
			}
		}
		else {
			m.insert(pair<string,fastq_read>(umi, fastq_read(seq, plus, qual)));
		}
	}
	in.close();
	return m;
}

int main(int argc, const char *argv[]) {
	if(argc != 5) {
		usage();
		exit(EXIT_FAILURE);
	}

	struct cmd_args args;
	args = parse_command_line(argc, argv);

	vector<string> fastq_files;
	fastq_files.push_back(args.ff);
	fastq_files.push_back(args.fr);

	map<string, fastq_read> mapper;
	for(int i = 0; i < fastq_files.size(); i++) {
		string curr_fq = fastq_files[i];
		mapper = process_fastq(curr_fq);
		write_to_file(mapper, base_name(curr_fq));
	}

	return 0;
}
