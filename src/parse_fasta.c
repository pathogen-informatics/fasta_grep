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

void search_for_query(char filename[], char ** search_queries, int number_of_queries)
{
	int l;
	int * found_queries;
	regex_t * regex_queries;
	
	regex_queries = (regex_t *)  malloc((number_of_queries+1)*sizeof(regex_t));
	found_queries  = (int *)  malloc((number_of_queries+1)*sizeof(int));
	int c= 0;
	for(c=0; c<number_of_queries;c++)
	{
		found_queries[c] = 0;
		regcomp(&regex_queries[c], search_queries[c], 0);
	}
	
	gzFile fp;
	kseq_t *seq; 
	
	fp = gzopen(filename, "r");
	seq = kseq_init(fp);
 
	int query_found = 0; 
	while ((l = kseq_read(seq)) >= 0) {
		int i = 0;
		for(i=0; i< number_of_queries; i++)
		{
			if(found_queries[i] == 1)
			{
				continue;
			}
			if(does_string_contain_query(seq->name.s, regex_queries[i]) == 1 )
			{
				printf(">%s\n", seq->name.s);
				printf("%s\n", seq->seq.s);
				query_found++;
				found_queries[i] = 1;
			}
			if(query_found == number_of_queries)
			{
				break;
			}
		}
		if(query_found == number_of_queries)
		{
			break;
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
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


