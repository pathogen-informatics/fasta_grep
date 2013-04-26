#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include "check_parse_fasta.h"
#include <regex.h>

START_TEST (check_does_string_contain_query)
{
	char *input_string = ">abc efg";
	char *input_query  = "(abc)";
	char *non_matching_input_query  = "zzz";
	
  regex_t  regex_input_query;
	regcomp(&regex_input_query, input_query, REG_EXTENDED);
	
	regex_t  regex_non_matching_input_query;
	regcomp(&regex_non_matching_input_query, non_matching_input_query, 0);
	
	fail_unless(does_string_contain_query(input_string,regex_input_query) == 1);
	fail_unless(does_string_contain_query(input_string,regex_non_matching_input_query) == 0);
	
}
END_TEST

START_TEST (check_parsing_fasta_files)
{
  char **input_queries;
	int i;
  input_queries = (char **) malloc((3+1)*sizeof(char *));
	for(i = 0; i < 3; i++)
	{
		input_queries[i] = (char *) malloc((40)*sizeof(char));
	}
	input_queries[0] = "non_matching_query";
  input_queries[1] = "matching_query";
  input_queries[2] = "another_matching";

	search_for_query("../tests/data/input_fasta_file", input_queries, 3);
	fail_unless(1==1);
}
END_TEST


Suite * parse_fasta_suite(void)
{
  Suite *s = suite_create("Parsing a fasta file");
  TCase *tc_parse_fasta = tcase_create("check_parsing_fasta_files");
  tcase_add_test(tc_parse_fasta, check_parsing_fasta_files);
	tcase_add_test(tc_parse_fasta, check_does_string_contain_query);
  suite_add_tcase(s, tc_parse_fasta);
  return s;
}
