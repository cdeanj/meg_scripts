#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <map>
#include "args.h"
#include "fastq_read.h"
using namespace std;

void write_to_file(const map<string,fastq_read> &mapper, const string &postfix) {
	ofstream out_match("match_" + postfix);
	ofstream out_mmatch("mismatch_" + postfix);
	map<string,fastq_read>::const_iterator it;
	for(it = mapper.begin(); it != mapper.end(); ++it) {
		if(!it->second._skip) {
			out_match << "@" << it->first << "/" << it->second._count << endl;
			out_match << it->second._seq << endl;
			out_match << it->second._plus << endl;
			out_match << it->second._qual << endl;
		}
		else {
			out_mmatch << "@" << it->first << endl;
			out_mmatch << it->second._seq << endl;
			out_mmatch << it->second._plus << endl;
			out_mmatch << it->second._qual << endl;
		}
	}
	out_match.close();
	out_mmatch.close();
}

string get_umi(const string &line) {
	int idx = line.find("|")+1;
	return line.substr(idx,24);
}

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
	cout << argc << endl;
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
		write_to_file(mapper, curr_fq);
	}

	return 0;
}
