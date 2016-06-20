#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <libgen.h>
#include <algorithm>
#include <assert.h>
#include <map>
#include <unordered_map>

#include "args.h"
#include "fastq_read.h"
using namespace std;

unordered_multimap<string,fastq_read> umm;

string base_name(string fp) {
	return string(basename(&fp[0]));
}

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

void transform(string &s1, const string &s2) {
	for(int i = 0; i < s1.length(); i++) {
		if(s1[i] != s2[i]) {
			s1[i] = 'N';
		}
	}
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
		int key_count = umm.count(umi);
		if(key_count >= 1) {
			if(seqs_are_equal(seq, umi)) {
				update_count(seq, umi);
			}
			else {
				umm.insert(make_pair(umi, fastq_read(seq, plus, qual)));
			}
		}
		else {
			umm.insert(make_pair(umi, fastq_read(seq, plus, qual)));
		}
	}
}

map<string,fastq_read> generate_consensus_fastq() {
	map<string,fastq_read> ordered_fastq;
	map<string,int> visited;
	for(auto parent_key = umm.begin(); parent_key != umm.end(); ++parent_key) {
		if(visited.count(parent_key->first) > 0) continue;
		int occ = 0;
		string umi = parent_key->first;
		string template_seq = parent_key->second._seq;
		auto key_range = umm.equal_range(parent_key->first);
		for(auto child_key = key_range.first; child_key != key_range.second; ++child_key) {
			visited[parent_key->first]++;
			if(template_seq != child_key->second._seq) {
				transform(template_seq, child_key->second._seq);
			}
			occ += child_key->second._count;
		}
		string plus = parent_key->second._plus;
		string qual = parent_key->second._qual;
		umi += ":" + to_string(occ);
		ordered_fastq.insert(make_pair(umi, fastq_read(template_seq, plus, qual)));
	}
	return ordered_fastq;
}

void write_fastq(const map<string,fastq_read> &mfq, const string &prefix, const string &basename) {
	ofstream ofs(prefix + "_" +  basename);
	for(auto it = mfq.begin(); it != mfq.end(); ++it) {
		ofs << it->first << '\n';
		ofs << it->second._seq << '\n';
		ofs << it->second._plus << '\n';
		ofs << it->second._qual << '\n';
	}
	ofs.close();
}


int main(int argc, const char *argv[]) {
	if(argc != 7) {
		usage();
		exit(EXIT_FAILURE);
	}
	struct cmd_args args;
	args = parse_command_line(argc, argv);

	vector<string> fastq_files;
	fastq_files.push_back(args.ff);
	fastq_files.push_back(args.fr);

	umm.reserve(50000000);
	map<string,fastq_read> mfq;
	for(int i = 0; i < fastq_files.size(); i++) {
		string curr_fq = fastq_files[i];
		process_fastq(curr_fq);
		mfq = generate_consensus_fastq();
		write_fastq(mfq, args.prefix, base_name(curr_fq));
		mfq.clear();
		umm.clear();
	}

	return 0;
}
