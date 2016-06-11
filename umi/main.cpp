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

struct match_and_mism {
	map<string,fastq_read> m;
	vector<pair<string,fastq_read> > mm;
};

string base_name(string fp) {
	return string(basename(&fp[0]));
}

struct match_and_mism get_m_and_mm(const multimap<string,fastq_read> &mm) {
	struct match_and_mism m_and_mm; 
	map<string,fastq_read> fm;
	vector<pair<string,fastq_read> > fmm;
	multimap<string,fastq_read>::const_iterator it;
	for(it = mm.begin(); it != mm.end(); ) {
		string k = it->first;
		if(mm.count(k) == 1) {
			fm.insert(pair<string,fastq_read>(k, it->second));
		}
		do {
			if(fm.count(k) == 0) {
				fmm.push_back(pair<string,fastq_read>(it->first, it->second));
			}
			++it;
		} while(it != mm.end() && k == it->first);
	}
	m_and_mm.m = fm;
	m_and_mm.mm = fmm;
	return m_and_mm;
}

void write_matches(const map<string,fastq_read> &m, const string &postfix) {
	ofstream ofs("matches_" + postfix);
	map<string,fastq_read>::const_iterator it;
	for(it = m.begin(); it != m.end(); ++it) {
		ofs << "@" << it->first << "/" << it->second._count << endl;
		ofs << it->second._seq << endl;
		ofs << it->second._plus << endl;
		ofs << it->second._qual << endl;
	}
	ofs.close();
}

void write_mismatches(const vector<pair<string,fastq_read> > &vp, const string &postfix) {
	ofstream ofs("mismatches_" + postfix);
	vector<pair<string,fastq_read> >::const_iterator it;
	for(it = vp.begin(); it != vp.end(); ++it) {
		ofs << "@" << it->first << "/" << it->second._count << endl;
		ofs << it->second._seq << endl;
		ofs << it->second._plus << endl;
		ofs << it->second._qual << endl;		
	}
	ofs.close();
}

void i_sect(map<string,fastq_read> &f, map<string,fastq_read> &r) {
	map<string,fastq_read>::iterator f_it, r_it;
	for(f_it = r.begin(); f_it != r.end(); ++f_it) {
		if(f.count(f_it->first) == 0) {
			r.erase(f_it->first);
		}
	}
	for(r_it = f.begin(); r_it != f.end(); ++r_it) {
		if(r.count(r_it->first) == 0) {
			f.erase(r_it->first);
		}
	}
}

void write_to_file(const vector<multimap<string, fastq_read> > &vmm, struct cmd_args args) {
	multimap<string,fastq_read> forward = vmm[0];
	multimap<string,fastq_read> reverse = vmm[1];
	struct match_and_mism f = get_m_and_mm(forward);
	struct match_and_mism r = get_m_and_mm(reverse);
	i_sect(f.m, r.m);
	string f_base_name = base_name(args.ff);
	string r_base_name = base_name(args.fr);
	write_matches(f.m, f_base_name);
	write_mismatches(f.mm, f_base_name);
	write_matches(r.m, r_base_name);
	write_mismatches(r.mm, r_base_name);
}

/*
 * Returns a 24-mer extracted from each each sequence identifier
*/
string get_umi(const string &line) {
	int idx = line.find("|")+1;
	return line.substr(idx,24);
}

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

	vector<multimap<string,fastq_read> > v;
	multimap<string, fastq_read> mapper;
	for(int i = 0; i < fastq_files.size(); i++) {
		string curr_fq = fastq_files[i];
		mapper = process_fastq(curr_fq);
		v.push_back(mapper);
	}
	write_to_file(v, args);

	return 0;
}
