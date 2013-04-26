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
#include <getopt.h>
#include <unistd.h>
#include "parse_fasta.h"
#include "../config.h"

const char* program_name;



void print_usage(FILE* stream, int exit_code)
{
  fprintf (stream, "Filter a FASTA file for based on a list of passed in search queries\n");
  fprintf (stream, "Usage:  %s [options] fasta_grep\n", program_name);
  fprintf (stream, "Version: %s\n", PACKAGE_VERSION);
  fprintf (stream,
           "  -f    fasta_file\n"
           "  -h    Display this usage information.\n\n"
);
  exit (exit_code);
}

int check_file_exists_or_exit(char * filename)
{
  if( access( filename, F_OK ) != -1 ) {
		return 1;
  } else {
		printf("Error: File '%s' doesnt exist\n",filename);
		print_usage(stderr, EXIT_FAILURE);
  }
}

int main (argc, argv) int argc; char **argv;
{
  int c;
  char fasta_filename[1024];

  program_name = argv[0];
  
  while (1)
    {
      static struct option long_options[] =
        {
					{"help",                no_argument,       0, 'h'},
          {"fasta_filename",      required_argument, 0, 'f'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;
      c = getopt_long (argc, argv, "hf:",
                       long_options, &option_index);
      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;
		    case 'h':
					print_usage(stdout, EXIT_SUCCESS);
	      case 'f':
	        strcpy(fasta_filename,optarg);
	        break;
        case '?':
          /* getopt_long already printed an error message. */
          break;
        default:
          print_usage(stdout, EXIT_SUCCESS);
          abort ();
        }
    }

		int i = 0;
		int number_of_searches =  argc - optind;
		
		if(number_of_searches == 0)
		{ 
    		print_usage(stdout, EXIT_FAILURE);	
		}
		
    char ** search_queries;
    search_queries = (char **) malloc((number_of_searches+1)*sizeof(char *));
		for(i = 0; i < number_of_searches; i++)
		{
			search_queries[i] = (char *) malloc((1024)*sizeof(char));
			strcpy(search_queries[i], argv[optind++]);
		}
		
		search_for_query(fasta_filename, search_queries, number_of_searches);
	  

    exit(EXIT_SUCCESS);
}
  


