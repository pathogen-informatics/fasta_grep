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


int number_of_valid_sequences(char filename[])
{
	int l;

	gzFile fp;
	kseq_t *seq; 
	
	fp = gzopen(filename, "r");
	seq = kseq_init(fp);
	
	int num_seqs = 0;
	
	while ((l = kseq_read(seq)) >= 0) {
		
	 if(does_sequence_have_start_or_stop_codons(seq->seq.s, seq->seq.l)  == 1)
	 {
		num_seqs++;
	 }
	}
	kseq_destroy(seq);
	gzclose(fp);
	
	return num_seqs;
}

void filter_out_invalid_sequences(char filename[])
{
	int l;

	int num_seqs = number_of_valid_sequences(filename);

	// initalise an to hold the names of the sequences and their lengths
	char ** sequence_names = (char **) malloc((num_seqs+1)*sizeof(char*));
	int * sequence_lengths = (int *) malloc((num_seqs+1)*sizeof(int));
	int j;
	for(j=0;j<num_seqs; j++)
	{
	  sequence_names[j] = (char *) malloc(1024*sizeof(char));
		sequence_lengths[j] = 0;
  }

	get_sequence_names_and_lengths( filename,num_seqs,sequence_names, sequence_lengths );
	flag_largest_sequence_if_duplicates( num_seqs,sequence_names, sequence_lengths );
	

	gzFile fp;
	kseq_t *seq; 
	
	fp = gzopen(filename, "r");
	seq = kseq_init(fp);
	int i =0;
	while ((l = kseq_read(seq)) >= 0) 
	{
	 if(does_sequence_have_start_or_stop_codons(seq->seq.s, seq->seq.l)  == 1)
	 {
		if(sequence_lengths[i] != -1)
		{
		  printf(">%s %s\n", seq->name.s, seq->comment.s);
	   	printf("%s\n", seq->seq.s);
  	}
		i++;
	 }
	}
	kseq_destroy(seq);
	gzclose(fp);
}

void flag_largest_sequence_if_duplicates(int num_seqs,char ** sequence_names,int * sequence_lengths )
{
	int i;
	char current_sequence_name[1024];
	strcpy(current_sequence_name, "");
	int current_largest_length = 0;
	int largest_index = 0;
	
	for(i = 0; i< num_seqs; i++)
	{
		if(strcmp(current_sequence_name, sequence_names[i]) == 0)
		{
			if(sequence_lengths[i] > current_largest_length)
			{
				sequence_lengths[largest_index] = -1;
				largest_index = i;
			}
			else
			{
				sequence_lengths[i] = -1;
			}
		}
		else
		{
			strcpy(current_sequence_name, sequence_names[i]);
			current_largest_length = sequence_lengths[i];
			largest_index = i;
		}
	}
}


void get_sequence_names_and_lengths(char filename[],int num_seqs,char ** sequence_names,int * sequence_lengths )
{
	int l;
	gzFile fp;
	kseq_t *seq; 
	
	fp = gzopen(filename, "r");
	seq = kseq_init(fp);

  // find the names and length of each sequence
	int i = 0;
	while ((l = kseq_read(seq)) >= 0) {
	 if(does_sequence_have_start_or_stop_codons(seq->seq.s, seq->seq.l)  == 1)
	 {
			strcpy(sequence_names[i],seq->name.s);
			sequence_lengths[i] = seq->seq.l;
			i++;
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

int does_sequence_have_start_or_stop_codons(char * input_string, int input_length)
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
	if(input_string[0] == 'M' || input_string[0] == 'V')
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


