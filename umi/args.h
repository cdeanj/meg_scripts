#ifndef _ARGS_H
#define _ARGS_H

#include <string>
#include <vector>

struct cmd_args {
	std::string ff;		// fastq forward reads file pointer
	std::string fr;		// fastq reverse reads file pointer
	std::string prefix;
};

void static usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: UMI \n");
	fprintf(stderr, "Contact: Chris Dean <cdean11@rams.colostate.edu>\n");
	fprintf(stderr, "./umi -first <forward> -second <reverse> -prefix <output prefix>\n");
	fprintf(stderr, "-first	STRING	fastq forward\n");
	fprintf(stderr, "-second STRING	fastq reverse\n\n");
}

static inline cmd_args
parse_command_line(int argc, const char *argv[]) {
	std::vector<std::string> arg_list(argv, argv+argc);
	cmd_args args;
	for(int i = 1; i < argc; i++) {
		if(arg_list[i] == "-first") 		args.ff = arg_list[++i];
		else if(arg_list[i] == "-second") 	args.fr = arg_list[++i];
		else if(arg_list[i] == "-prefix") 	args.prefix = arg_list[++i];
		else {
			usage();
			exit(EXIT_FAILURE);
		}
	}
	return args;
}

#endif // _ARGS_H
