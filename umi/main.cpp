#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <libgen.h>
#include <string.h>
#include <algorithm>
#include <chrono>
#include <map>
#include <unordered_map>

#include "args.h"
#include "fastq_read.h"
using namespace std;

unordered_multimap<string,fastq_read> umm;

struct match_and_mism {
	unordered_map<string,fastq_read> u_m;
	vector<pair<string,fastq_read> > vp_mm;
};

string base_name(string fp) {
	return string(basename(&fp[0]));
}

struct match_and_mism get_m_and_mm(const unordered_multimap<string,fastq_read> &umm_reads) {
	struct match_and_mism m_and_mm; 
	unordered_map<string,fastq_read> fm;
	vector<pair<string,fastq_read> > fmm;
	unordered_multimap<string,fastq_read>::const_iterator it;
	for(it = umm_reads.begin(); it != umm_reads.end(); ) {
		string k = it->first;
		if(umm_reads.count(k) == 1) {
			fm.insert(pair<string,fastq_read>(k, it->second));
		}
		do {
			if(fm.count(k) == 0) {
				fmm.push_back(pair<string,fastq_read>(it->first, it->second));
			}
			++it;
		} while(it != umm_reads.end() && k == it->first);
	}
	m_and_mm.u_m = fm;
	m_and_mm.vp_mm = fmm;
	return m_and_mm;
}

void write_matches(const map<string,fastq_read> &m, const string &postfix) {
	ofstream ofs("matches_" + postfix);
	for(auto it = m.begin(); it != m.end(); ++it) {
		ofs << "@" << it->first << "/" << it->second._count << '\n';
		ofs << it->second._seq << '\n';
		ofs << it->second._plus << '\n';
		ofs << it->second._qual << '\n';
	}
	ofs.close();
}

void write_mismatches(const vector<pair<string,fastq_read> > &vp, const string &postfix) {
	ofstream ofs("mismatches_" + postfix);
	for(auto it = vp.begin(); it != vp.end(); ++it) {
		ofs << "@" << it->first << "/" << it->second._count << '\n';
		ofs << it->second._seq << '\n';
		ofs << it->second._plus << '\n';
		ofs << it->second._qual << '\n';		
	}
	ofs.close();
}

void i_sect(unordered_map<string,fastq_read> &f, unordered_map<string,fastq_read> &r) {
	auto f_it = r.begin();
	while(f_it != r.end()) {
		if(f.count(f_it->first) == 0) {
			f_it = r.erase(f_it);
		}
		else {
			++f_it;
		}
	} 
	auto r_it = f.begin();
	while(r_it != f.end()) {
		if(r.count(r_it->first) == 0) {
			r_it = f.erase(r_it);
		}
		else {
			++r_it;
		}
	}	
}

void write_to_file(const vector<unordered_multimap<string, fastq_read> > &vmm, struct cmd_args args) {
	struct match_and_mism f = get_m_and_mm(vmm[0]);
	struct match_and_mism r = get_m_and_mm(vmm[1]);
	i_sect(f.u_m, r.u_m);
	map<string,fastq_read> s_f(f.u_m.begin(), f.u_m.end());
	map<string,fastq_read> s_r(r.u_m.begin(), r.u_m.end());

	string f_base_name = base_name(args.ff);
	string r_base_name = base_name(args.fr);

	write_matches(s_f, f_base_name);
	write_mismatches(f.vp_mm, f_base_name);
	write_matches(s_r, r_base_name);
	write_mismatches(r.vp_mm, r_base_name);
}

/*
 * Returns a 24-mer extracted from each each sequence identifier
*/
string get_umi(const string &line) {
	int idx = line.find("|")+1;
	return line.substr(idx,24);
}

bool seqs_are_equal(const string &seq, const string &umi) {
	auto key_range = umm.equal_range(umi);
	for(auto it = key_range.first; it != key_range.second; ++it) {
		if(it->second._seq == seq) {
			return true;
		}
	}
	return false;
}

void update_count(const string &seq, const string &umi) {
	auto key_range = umm.equal_range(umi);
	for(auto it = key_range.first; it != key_range.second; ++it) {
		if(it->second._seq == seq) {
			it->second._count += 1;
			break;
		}
	}
}

bool seq_was_added(const string &seq, const string &umi) {
	auto r = umm.equal_range(umi);
	for(auto d = r.first; d != r.second; ++d) {
		if(d->second._seq == seq) {
			return true;
		}
	}
	return false;
}

void process_fastq(const string &f) {
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
		int c = umm.count(umi);
		if(c == 1) {
			// is the current sequence identical to the sequence associated with the umi?
			if(seqs_are_equal(seq, umi)) {
				update_count(seq, umi);	
			}
			// if these sequences are not unique, then we should add this sequence to the current list
			else {
				umm.insert(pair<string, fastq_read>(umi, fastq_read(seq, plus, qual)));
			}
		}
		// if the amount of keys equal to umi is greater than one, then there are more than one sequences corresponding to this umi
		// as a result, we update the count of the sequence associated to that some umi
		else if(c > 1) {
			if(seq_was_added(seq, umi)) {
				update_count(seq, umi);
			}
			else {
				umm.insert(pair<string, fastq_read>(umi, fastq_read(seq, plus, qual)));
			}
		}
		// umi doesn't exist yet, let's add it
		else {
			umm.insert(pair<string, fastq_read>(umi, fastq_read(seq, plus, qual)));
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

	umm.reserve(100000000);

	vector<unordered_multimap<string,fastq_read> > vmm;
	multimap<string, fastq_read> mapper;
	for(int i = 0; i < fastq_files.size(); i++) {
		string curr_fq = fastq_files[i];
		process_fastq(curr_fq);
		vmm.push_back(umm);
		umm.clear();
	}
	write_to_file(vmm, args);

	return 0;
}
