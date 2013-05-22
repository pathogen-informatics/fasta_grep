/*
 *  Wellcome Trust Sanger Institute
 *  Copyright (C) 2013  Wellcome Trust Sanger Institute
 *  
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <regex.h>
#include <sys/types.h>
#include "kseq.h"
#include "parse_fasta.h"

KSEQ_INIT(gzFile, gzread)

void filter_out_invalid_sequences(char filename[])
{
	int l;
	
	regex_t  regex_query;

	gzFile fp;
	kseq_t *seq; 
	
	fp = gzopen(filename, "r");
	seq = kseq_init(fp);
 
	while ((l = kseq_read(seq)) >= 0) {
		int i = 0;
		
	 if(does_sequence_have_start_and_stop_codons(seq->seq.s, seq->seq.l)  == 1)
	 {
	 	printf(">%s %s\n", seq->name.s, seq->comment.s);
	 	printf("%s\n", seq->seq.s);
	 }
	}
	kseq_destroy(seq);
	gzclose(fp);
}





void search_for_query(char filename[], char ** search_queries, int number_of_queries)
{
	int l;
	int * found_queries;
	
	regex_t  regex_query;

	
	// Compile a single regular expression once
	char * combined_search_query;
	combined_search_query =  (char *) malloc(((number_of_queries*1024) +1)*sizeof(char));
	strcpy(combined_search_query,"");
	strcat(combined_search_query,"(");
	int c =0;
	for(c = 0; c< number_of_queries; c++)
	{
		if(c > 0)
		{
			strcat(combined_search_query,"|");
		}
		strcat(combined_search_query,search_queries[c]);
	}
	strcat(combined_search_query,")");
	regcomp(&regex_query, combined_search_query, REG_EXTENDED);
	
	gzFile fp;
	kseq_t *seq; 
	
	fp = gzopen(filename, "r");
	seq = kseq_init(fp);
 
	int query_found = 0; 
	while ((l = kseq_read(seq)) >= 0) {
		int i = 0;
		
	 if(does_string_contain_query(seq->name.s, regex_query) == 1)
	 {
	 	printf(">%s %s\n", seq->name.s, seq->comment.s);
	 	printf("%s\n", seq->seq.s);
	 	query_found++;
	 }
	 if(query_found == number_of_queries)
	 {
	 	break;
	 }

	}
	kseq_destroy(seq);
	gzclose(fp);
}

int does_sequence_have_start_and_stop_codons(char * input_string, int input_length)
{
	// There shouldnt be a stop codon in the middle of the sequence
	int i;
	for(i = 0; i < input_length-1; i++)
	{
		if(input_string[i] == '*')
		{
			return 0;
		}
	}
	
	// Should usually be a start codon
	int start_codon_found = 0;
	if(input_string[0] == 'M')
	{
		start_codon_found = 1;
	}
	
	
	// Should usually be an end codon
	int end_codon_found = 0;
	if(input_string[input_length-1] == '*')
	{
		end_codon_found = 1;
	}
	
	if(start_codon_found == 1 || end_codon_found == 1)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int does_string_contain_query(char * input_string, regex_t regex)
{
	int reti;
	/* Execute regular expression */
	reti = regexec(&regex, input_string, 0, NULL, 0);

	if( !reti ){
		return 1;
	}
	else if( reti == REG_NOMATCH ){
		return 0;
	}	
	return 0;
}


