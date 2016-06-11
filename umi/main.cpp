#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <libgen.h>
#include <string.h>
#include <algorithm>
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
void write_to_file(const multimap<string,fastq_read> &mm, const string &postfix) {
	ofstream ofs_m("match_" + postfix);
	ofstream ofs_mm("mismatch_" + postfix);

	for(auto it = mm.begin(); it != mm.end(); ++it) {
		if(mm.count(it->first) == 1) {
			ofs_m << "@" << it->first << "/" << it->second._count << endl;
			ofs_m << it->second._seq << endl;
			ofs_m << it->second._plus << endl;
			ofs_m << it->second._qual << endl;		
		}
		else {
			ofs_mm << "@" << it->first << "/" << it->second._count << endl;
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
/*map<string, fastq_read> process_fastq(const string &f) {
	map<string, vector<fastq_read> > m;
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
}*/

bool seqs_are_equal(const multimap<string, fastq_read> &mm, const string &seq, const string &umi) {
	auto key_range = mm.equal_range(umi);
	for(auto it = key_range.first; it != key_range.second; ++it) {
		if(it->second._seq == seq) {
			return true;
		}
	}
	return false;
}

void update_count(multimap<string, fastq_read> &mm, const string &seq, const string &umi) {
	auto key_range = mm.equal_range(umi);
	for(auto it = key_range.first; it != key_range.second; ++it) {
		if(it->second._seq == seq) {
			it->second._count += 1;
			break;
		}
	}
}

bool seq_was_added(const multimap<string, fastq_read> mm, const string &seq, const string &umi, const string &qual) {
	fastq_read fr(seq, "+", qual);
	auto const &r = mm.equal_range(umi);
	bool found = std::any_of(r.first, r.second,
				[&fr](decltype(mm)::value_type const &p) {
					return p.second._seq == fr._seq;
				});
	return found;
}

multimap<string, fastq_read> process_fastq(const string &f) {
	multimap<string, fastq_read> mm;
	string line, umi, seq, plus, qual;

	ifstream in(f);
	if(!in) {
		cerr << "Could not open fasta file" << endl;
		exit(EXIT_FAILURE);
	}
	while(getline(in, line)) {
		if(line[0] == '@') {
			umi = get_umi(line);
			getline(in, seq);
			getline(in, plus);
			getline(in, qual);
		}
		// if the umi count is one, then we know we're dealing with umi's with identical sequences (so far)
		if(mm.count(umi) == 1) {
			// is the current sequence identical to the sequence associated with the umi?
			if(seqs_are_equal(mm, seq, umi)) {
				update_count(mm, seq, umi);	
			}
			// if these sequences are not unique, then we should add this sequence to the current list
			else {
				mm.insert(pair<string, fastq_read>(umi, fastq_read(seq, plus, qual)));
			}
		}
		// if the amount of keys equal to umi is greater than one, then there are more than one sequences corresponding to this umi
		// as a result, we update the count of the sequence associated to that some umi
		else if(mm.count(umi) > 1) {
			if(seq_was_added(mm, seq, umi, qual)) {
				update_count(mm, seq, umi);
			}
			else {
				mm.insert(pair<string, fastq_read>(umi, fastq_read(seq, plus, qual)));
			}
		}
		// umi doesn't exist yet, let's add it
		else {
			mm.insert(pair<string, fastq_read>(umi, fastq_read(seq, plus, qual)));
		}
	}
	return mm;	
}

void print_mm(multimap<string, fastq_read> &mm) {
	for(auto it = mm.begin(); it != mm.end(); ++it) {
		if(mm.count(it->first) == 1) {
			cout << "MATCH" << endl;
			cout << "@" << it->first << "/" << it->second._count << endl;
			cout << it->second._seq << endl;
			cout << it->second._plus << endl;
			cout << it->second._qual << endl;
		}
		else {
			cout << "MISMATCH" << endl;
			cout << "@" << it->first << "/" << it->second._count << endl;
			cout << it->second._seq << endl;
			cout << it->second._plus << endl;
			cout << it->second._qual << endl;
		}
	}
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

	multimap<string, fastq_read> mapper;
	for(int i = 0; i < fastq_files.size(); i++) {
		string curr_fq = fastq_files[i];
		mapper = process_fastq(curr_fq);
		//print_mm(mapper);	
		mapper = process_fastq(curr_fq);
		write_to_file(mapper, base_name(curr_fq));
	}

	return 0;
}
