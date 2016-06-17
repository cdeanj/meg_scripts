#include <iostream>
#include <fstream>
#include <string>
#include <string>
#include <sstream>
#include <libgen.h>
#include <map>
using namespace std;

int s_to_i(const string &s) {
	istringstream ss(s);
	int i;
	ss >> i;
	return i;
}

map<int, int> get_umi_count(const string &fp) {
	ifstream ifs(fp);
	if(!ifs) {
		cerr << "Could not open fastq file" << endl;
		exit(EXIT_FAILURE);
	}
	int member, colon_idx;
	map<int, int> mapper;
	string umi, seq, plus, qual, occ;
	while(getline(ifs, umi)) {
		colon_idx = umi.find(":");
		occ = umi.substr(colon_idx+1, string::npos);
		member = s_to_i(occ);
		mapper[member]++;
		getline(ifs, seq);
		getline(ifs, plus);
		getline(ifs, qual);
	}
	ifs.close();
	return mapper;
}

string base_name(string fp) {
	string b_name = string(basename(&fp[0]));
	int ext_idx = b_name.find_last_of(".");
	return b_name.substr(0, ext_idx);
}

void write_umi_count(const map<int, int> &m, const string &prefix, const string &basename) {
	ofstream ofs(prefix + "_" + basename);
	map<int, int>::const_iterator it;
	for(it = m.begin(); it != m.end(); ++it) {
		ofs << it->first << "," << it->second << '\n';
	}
	ofs.close();
}

int main(int argc, char *argv[]) {
	if(argc != 4) {
		cerr << "Usage: <forward> <reverse> <output_prefix>" << endl;
		return -1;
	}

	string forward_fp = argv[1];
	string reverse_fp = argv[2];
	string output_prefix = argv[3];

	map<int, int> forward_count = get_umi_count(forward_fp);
	map<int, int> reverse_count = get_umi_count(reverse_fp);

	string forward_base_name = base_name(forward_fp);
	string reverse_base_name = base_name(reverse_fp);

	write_umi_count(forward_count, output_prefix, forward_base_name);
	write_umi_count(reverse_count, output_prefix, reverse_base_name);

	return 0;
}
