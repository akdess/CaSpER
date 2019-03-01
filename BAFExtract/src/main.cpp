#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "ansi_string.h"
#include "utils.h"
#include "mapped_read_tools.h"
#include "genomics_coords.h"
//#include "genome_sequence_tools.h"
#include <algorithm>

using namespace std;

bool __DUMP_MAPPED_READ_TOOLS_MSGS__ = false;

unsigned short** allocate_pileup(int l_sig)
{
  unsigned short** loaded_pileup = new unsigned short*[5];
  for (int i_allele = 0; i_allele < 5; i_allele++)
  {
    loaded_pileup[i_allele] = new unsigned short[l_sig + 2];
    memset(loaded_pileup[i_allele], 0, sizeof(unsigned short) * (l_sig + 2));
  } // i_allele loop.

  return loaded_pileup;
}
int* get_per_char_as_nuc_2_num_coding_array()
{
  // code all non-nucleotide values as 4, others as a:0, c:1, g:2, t:3.
  int* per_nuc_val = new int[256];
  for(int i = 0; i < 256; i++)
  {
    per_nuc_val[i] = 4;
  } // i loop.
  per_nuc_val[(int)'A'] = 0;
  per_nuc_val[(int)'C'] = 1;
  per_nuc_val[(int)'G'] = 2;
  per_nuc_val[(int)'T'] = 3;
  per_nuc_val[(int)'U'] = 3;
  per_nuc_val[(int)'a'] = 0;
  per_nuc_val[(int)'c'] = 1;
  per_nuc_val[(int)'g'] = 2;
  per_nuc_val[(int)'t'] = 3;
  per_nuc_val[(int)'u'] = 3;

  return(per_nuc_val);
}

void normalize_chr_id(char* raw_chr_id)
{
  if(t_string::starts_with(raw_chr_id, "chr") || 
    t_string::starts_with(raw_chr_id, "Chr") ||
    t_string::starts_with(raw_chr_id, "CHR"))
  {
    // Get rid of the first 3 characters.
    char temp[1000];
    strcpy(temp, &raw_chr_id[3]);
    strcpy(raw_chr_id, temp);
  }
  else
  {
    // No need to change.
    return;
  }
}

void load_chromosome_lengths_per_tabbed_file(char* chr_lengths_fp, vector<char*>* chr_ids, vector<int>* chr_lengths)
{
  FILE* f_l_chrs = open_f(chr_lengths_fp, "r");

  while(1)
  {
    char* cur_line = getline(f_l_chrs);
    if(cur_line == NULL)
    {
      break;
    }

    char cur_chr_id[1000];
    int cur_l_chr = 0;
    if(sscanf(cur_line, "%s %d", cur_chr_id, &cur_l_chr) != 2)
    {
      fprintf(stderr, "Could not read chromosome lengths from: %s\n", cur_line);
      exit(0);
    } // check if reading is succesful.

    // Normalize the chromosome id before adding it to the list.
    normalize_chr_id(cur_chr_id);
    chr_ids->push_back(t_string::copy_me_str(cur_chr_id));
    chr_lengths->push_back(cur_l_chr);
  } // file reading loop.

  fclose(f_l_chrs);
}



unsigned short** load_compressed_pileups(char* cur_comp_allele_fp, int& l_pileup)
{
  FILE* f_comp = NULL;
  if(t_string::ends_with(cur_comp_allele_fp, ".bin"))
  {
    f_comp = open_f(cur_comp_allele_fp, "rb");
  }
  else if(t_string::ends_with(cur_comp_allele_fp, ".bin.gz"))
  {
    char ungzip_cmd[1000];
    sprintf(ungzip_cmd, "gzip -cd %s", cur_comp_allele_fp);
#ifdef WIN32
    f_comp = _popen(ungzip_cmd, "r");
#else 
    f_comp = popen(ungzip_cmd, "r");
#endif
  }

  if(f_comp == NULL)
  {
    fprintf(stderr, "Could not open %s\n", cur_comp_allele_fp);
    exit(0);
  }

  int l_read_sig = 0;
  fread(&l_read_sig, sizeof(int), 1, f_comp);
  int l_sig = l_read_sig;

  l_pileup = l_sig;

  fprintf(stderr, "Loading pileup of length %d\n", l_sig);

  unsigned short** loaded_pileup = allocate_pileup(l_sig);

  // Go over all the positions.
  int i = 1;
  while(i <= l_sig)
  {
    // Read the existence flag.
    unsigned char existence_flag = 0;
    fread(&existence_flag, sizeof(unsigned char), 1, f_comp);


    // Check an RLE case.
    if(existence_flag == 0xFF)
    {
      unsigned int l_RLE = 0;
      fread(&l_RLE, sizeof(unsigned int), 1, f_comp);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
      fprintf(stderr, "Loading RLE of length %d @ %d\n", l_RLE, i);

      // When we add this, we move upto the location where a 0-run ends.
      i += l_RLE;
    } // RLE check.
    else
    {
      // Read the first byte, parse the data range and which alleles exist.
      unsigned char drange_flag = ((existence_flag & (1 << 5)) >> 5);

      for(int allele_i = 0; allele_i < 5; allele_i++)
      {
        // Do we have an entry here?
        if((existence_flag & (1 << allele_i)) > 0)
        {
          if(drange_flag == 1)
          { 
            unsigned short current_count_short = 0;
            fread(&(current_count_short), sizeof(unsigned short), 1, f_comp);
            loaded_pileup[allele_i][i] = current_count_short;
          }
          else
          {
            unsigned char current_count_char = 0;
            fread(&(current_count_char), sizeof(unsigned char), 1, f_comp);
            loaded_pileup[allele_i][i] = (unsigned short)current_count_char;
          }
        } // count check for current position.
      } // allele_i loop.

      i++;
    }
  } // i loop.

  if(t_string::ends_with(cur_comp_allele_fp, ".bin"))
  {
    fclose(f_comp);
  }
  else if(t_string::ends_with(cur_comp_allele_fp, ".bin.gz"))
  {
#ifdef WIN32
    _pclose(f_comp);
#else 
    pclose(f_comp);
#endif
  }

  // Free memory.
  return(loaded_pileup);
}


enum{VAL, TYPE};
bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced)
{
  int i = 0;

  is_read_spliced = false;

  int state = VAL;
  while(mapping_map_str[i] != 0)
  {   
    if(state == VAL)
    {
      //fprintf(stderr, "%c (%d)\n", quality_str[i], 0);
      // MIDNSHPX=
      if(mapping_map_str[i] == 'M' ||
        mapping_map_str[i] == 'I' ||
        mapping_map_str[i] == 'D' ||
        mapping_map_str[i] == 'N' ||
        mapping_map_str[i] == 'S' ||    
        mapping_map_str[i] == 'H' ||
        mapping_map_str[i] == 'P' ||
        mapping_map_str[i] == 'X' ||
        mapping_map_str[i] == '=')
      {
        state = TYPE;

        if(mapping_map_str[i] != 'M')
        {
          is_read_spliced = true;
        }
      }
      else if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
      {
        // State is still VAL.
      }
      else
      {
        return(false);
      }
    }
    else if(state == TYPE)
    {
      //fprintf(stderr, "%c (%d)\n", quality_str[i], 1);
      // A number is expected.
      if(mapping_map_str[i] >= '0' && mapping_map_str[i] <= '9')
      {
        state = VAL;
      }
      else
      {
        return(false);
      }
    }

    // Move to next character.
    i++;
  }

  return(true);
}


void get_next_entry_per_mapp_map_string(char* mapping_map_str,
                    int& i_mapp_map, 
                    bool& is_matching,
                    //t_string* cur_entry_length_str,
                    int& l_cur_entry,
                    char& entry_type_char)
{ 
  // Clean the length string.
  //cur_entry_length_str->empty();
  l_cur_entry = 0;

  // Get the next entry in the cigar string.
  while(mapping_map_str[i_mapp_map] != 0)
  {
    if(mapping_map_str[i_mapp_map] < '0' || mapping_map_str[i_mapp_map] > '9')
    {
      break;
    }
    //cur_entry_length_str->concat_char(mapping_map_str[i_mapp_map]);
    l_cur_entry = l_cur_entry*10 + (int)(mapping_map_str[i_mapp_map]-'0');
    i_mapp_map++;
  }

  is_matching = false;
  if(mapping_map_str[i_mapp_map] == 'M')
  {
    //fprintf(stderr, "Adding matching length of %d\n", l_cur_entry);
    is_matching = true;
  }
  else
  {
    //fprintf(stderr, "Adding some other length of %d\n", l_cur_entry);
  } 

  entry_type_char = mapping_map_str[i_mapp_map];

  // Move over the current entry identifier.
  i_mapp_map++;
}

void compress_nucleotide_pileup_track(unsigned short** pileup, int l_sig, char* op_fp)
{
  FILE* f_op = open_f(op_fp, "wb");
  
  // Write the length of pileup
  fwrite(&l_sig, sizeof(int), 1, f_op);
  int n_4bit_posns = 0;
  int n_8bit_posns = 0;
  int n_12bit_posns = 0;
  int n_14bit_posns = 0;
  int n_16bit_posns = 0;
  int n_0_signal = 0;

  // This is the 0 RLE start position.
  int cur_0_run_start = -1;
  for(int i = 1; i <= l_sig; i++)
  {
    // Dump the current 5 levels:     
    // Generate the existence flag.
    //fprintf(stderr, "Position %d:\n", i);
    unsigned char existence_flag = 0;
    unsigned char drange_flag = 0;
    unsigned short max_val = 0;
    unsigned int total_val = 0;
    for(int allele_i = 0; allele_i < 5; allele_i++)
    {
      if(pileup[allele_i][i] > 0)
      {
        total_val += (unsigned int)(pileup[allele_i][i]);

        if(max_val < pileup[allele_i][i])
        {
          max_val = pileup[allele_i][i];
        }

        existence_flag = existence_flag | (1 << allele_i);
        //fprintf(stderr, "Nuc %d: %d\n", allele_i, pileup[allele_i][i]);
        //getc(stdin);
      } // signal check.
    } // allele_i loop.

    // Write the 1 byte phred-based partitioning flag.
    // Following is an RLE check.
    if(total_val == 0)
    {
      if(cur_0_run_start == -1)
      {
        // Initiate the 0 run.
        cur_0_run_start = i;

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
        fprintf(stderr, "Initiated a 0-run @ %d\n", i);
      }
      else
      {
        // We are already in a zero run.
      }
    } // total_val 0 check.
    else
    {
      // Check if we will do a RLE dump.
      // Make sure there was a run of zeros before this position.
      if(cur_0_run_start != -1)
      {
        // Check the RLE length and dump if requested.
        bool do_RLE_dump = ((i - cur_0_run_start) > 5);

        if(do_RLE_dump)
        {
          //fprintf(stderr, "RLE dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

          // Do an RLE dump: Note that an RLE dump is guaranteed to be more efficient as long as the length run is longer than 5.
          // Note that the top 2 bits are reserved for RLE, setting the RLE indicator to FF.
          unsigned char RLE_indicator = 0xFF;
          fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
          fprintf(stderr, "Writing: %x\n", RLE_indicator);

          // Write the length of 0 length region: Note that i does not have a 0, it must not be included.
          unsigned int l_RLE = (i - cur_0_run_start);
          fwrite(&l_RLE, sizeof(unsigned int), 1, f_op);

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
          fprintf(stderr, "Writing l_RLE: %d\n", l_RLE);
        } // RLE dump check.
        else
        {
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
          fprintf(stderr, "Per element dumping a 0-run @ %d till %d\n", cur_0_run_start, i);

          // Save the position with normal dump where we dump 0's.
          for(int j = cur_0_run_start; j < i; j++)
          {
            // Write 0 to each position.
            unsigned char RLE_indicator = 0;
            fwrite(&RLE_indicator, sizeof(unsigned char), 1, f_op);
          } // j loop.
        } // non-RLE check.

        // Reset the RLE start position.
        cur_0_run_start = -1;
      } // 0-run check. 

      if(max_val == 0)
      {
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
        fprintf(stderr, "Sanity check failed: Max value is 0 on a non-zero dumping path.\n");
        exit(0);
      }
      
      // Set the data range flag: 1 for short, 0 for char.
      if(max_val >= (1 << 8)-1)
      {
        drange_flag = 1;
      }

      // If val is not char, flag it to make sure we dump the right number of bytes.
      existence_flag = existence_flag | (drange_flag << 5);

      // Write existence flag.
      fwrite(&existence_flag, sizeof(unsigned char), 1, f_op);

      // Dump the numbers.
      for(int allele_i = 0; allele_i < 5; allele_i++)
      {
        if(pileup[allele_i][i] > 0)
        {
          if(drange_flag == 1)
          { 
            unsigned short current_count_short = (unsigned short)(pileup[allele_i][i]);
            fwrite(&(current_count_short), sizeof(unsigned short), 1, f_op);
          }
          else
          {
            unsigned char current_count_char = (unsigned char)(pileup[allele_i][i]);
            fwrite(&(current_count_char), sizeof(unsigned char), 1, f_op);
          }
        } // count check for current position.
      } // allele_i loop.

    } // total value check.

    // Get some simple stats.
    if(max_val == 0)
    {
      n_0_signal++;
    }
    else if(max_val <= ((2<<4)-1))
    {
      n_4bit_posns++;
    }
    else if(max_val <= ((2<<8)-1))
    {
      n_8bit_posns++;
    }
    else if(max_val <= ((2<<12)-1))
    {
      n_12bit_posns++;
    }
    else if(max_val <= ((2<<14)-1))
    {
      n_14bit_posns++;
    }
    else if(max_val <= ((2<<16)-1))
    {
      n_16bit_posns++;
    }
  } // i loop.
  fclose(f_op);
  fprintf(stderr, "%d 0 signal, %d 4-bit, %d 8-bit, %d 12-bit, %d 14-bit, %d 16-bit positions.\n", n_0_signal, n_4bit_posns, n_8bit_posns, n_12bit_posns, n_14bit_posns, n_16bit_posns);

  // Now, load and test if we received the correct values.
  bool loading_check = false;
  if(loading_check)
  {
    fprintf(stderr, "Loading and comparing.\n");

    // Load the compressed pileup; then do sanity check.
    int l_loaded_pileup = 0;
    unsigned short** loaded_pileup = load_compressed_pileups(op_fp, l_loaded_pileup);

    // Compare the loaded and the actual pileups.
    for(int i_allele = 0; i_allele < 5; i_allele++)
    {
      for(int i = 1; i <= l_sig; i++)
      {
        if(loaded_pileup[i_allele][i] != pileup[i_allele][i])
        {
          fprintf(stderr, "Sanity check failed: posn: %d; allele: %d: %d,%d\n", 
            i, i_allele, loaded_pileup[i_allele][i], pileup[i_allele][i]);
          exit(0);
        }
      } // i loop.        
    } // i_allele loop.

    // Free memory for both loaded and dumped pileups.
    for(int i_allele = 0; i_allele < 5; i_allele++)
    {
      delete [] loaded_pileup[i_allele];
    } // i_allele loop.
  }
  else
  {
    fprintf(stderr, "Skipping loading check.\n");
  }
}


//TODO: Following may not be complete up to the SAM file formatting's cigar string characters, must be checked.
bool check_genome_index_update_per_CIGAR_entry(char entry_char)
{
	//if(entry_char == 'D' || 
	//	entry_char == 'M' ||
	//	entry_char == 'N' ||
	//	entry_char == 'H')
	if (entry_char == 'D' ||
		entry_char == 'M' ||
		entry_char == 'N')
	{
		return(true);
	}

	return(false);
}

//TODO: Following may not be complete up to the SAM file formatting's cigar string characters, must be checked.
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char)
{
	//if(entry_char == 'H' ||
	//if(entry_char == 'P' ||
	//if(entry_char == 'X' ||
	if (entry_char == 'S' ||
		entry_char == 'M' ||
		entry_char == 'I')
	{
		return(true);
	}

	return(false);
}


void dump_nucleotide_pileup_per_SAM_file(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, int min_phred_qual, int& n_processed_reads)
{
  fprintf(stderr, "Generating pileup from SAM file with min mapp qual %d, min base qual %d\n", min_mapp_qual, min_phred_qual);

  //// Init the file.
  //char summary_fp[1000];
  //sprintf(summary_fp, "%s/%s", op_dir, t_config_params::OP_filenames[OP_PILEUP_SUMMARY_FN]);
  //FILE* f_summary = open_f(summary_fp, "w");
  //fclose(f_summary);

  // Initialize the number of processed readss.
  n_processed_reads = 0;

  //// Make sure we do not overwrite on an existing set of pileup file, since they are time consuming to generate.
  //bool allele_counts_are_there = true;
  //for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
  //{
  //  char cur_allele_fp[1000];
  //  sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
  //  if(!check_file(cur_allele_fp))
  //  {
  //    // This file does not exist, not all the allele counts are there.
  //    allele_counts_are_there = false;
  //    break;
  //  }
  //} // i_chr loop.

  //if(allele_counts_are_there)
  //{
  //  fprintf(stderr, "All the chromosome pileups exist, will not process data.\n");
  //  return;
  //}

  fprintf(stderr, "Dumping the pileups per SAM file %s\n", sam_fp);

  fprintf(stderr, "Allocating pileup memory.\n");
  unsigned short*** per_chrom_nuc_count_per_allele = new unsigned short**[(int)chr_ids->size()];
  //int** per_chrom_coverage = new int*[(int)chr_ids->size()];
  for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
  {
    int l_sig = chr_lengths->at(i_chr);

    // Check the existing pileup file, if it is there, add to it.
    char cur_allele_fp[1000];
    sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));
    if(check_file(cur_allele_fp))
    {
      int l_loaded_pileup = 0;
      per_chrom_nuc_count_per_allele[i_chr] = load_compressed_pileups(cur_allele_fp, l_loaded_pileup);
      if(l_loaded_pileup != chr_lengths->at(i_chr))
      {
        fprintf(stderr, "Loaded pileup is not the same length as the chromosome info: %s: %d, %d\n", chr_ids->at(i_chr), chr_lengths->at(i_chr));
        exit(0);
      }
      else
      {
        fprintf(stderr, "Loaded existing pileup from %s.\n", cur_allele_fp);
      }
    }
    else
    {
      per_chrom_nuc_count_per_allele[i_chr] = allocate_pileup(l_sig);
    }
  } // i_chr loop.
  fprintf(stderr, "Done.\n");

  int max_n_alleles_per_posn = 1<<16;
  fprintf(stderr, "Maximum # of alleles at each position: %d\n", max_n_alleles_per_posn);

  // Start reading the file.
  FILE* f_sam = NULL;
  if(strcmp(sam_fp, "stdin") == 0)
  {
    f_sam = stdin;
  }
  else
  {
    f_sam = fopen(sam_fp, "r");
  }

  int* per_nuc_val = get_per_char_as_nuc_2_num_coding_array();
  //for(int i = 0; i < 256; i++)
  //{
  //  per_nuc_val[i] = 4;
  //} // i loop.
  //per_nuc_val[(int)'A'] = 0;
  //per_nuc_val[(int)'C'] = 1;
  //per_nuc_val[(int)'G'] = 2;
  //per_nuc_val[(int)'T'] = 3;
  //per_nuc_val[(int)'U'] = 3;
  //per_nuc_val[(int)'a'] = 0;
  //per_nuc_val[(int)'c'] = 1;
  //per_nuc_val[(int)'g'] = 2;
  //per_nuc_val[(int)'t'] = 3;
  //per_nuc_val[(int)'u'] = 3;

  // Enter file reading loop.
  char read_id[1000];
  //int flag;
  char chrom[100];
  int chr_index;
  char mapping_map_str[10000];
  char cur_fragment[1000];
  char flag_str[100];
  char _chr_index_str[100];
  char phred_quality_str[100000];
  char mapp_quality_str[100];

  int n_unmapped_reads = 0;
  int n_low_quality_reads = 0;

  int phred_score_base = -123;

  while(1)
  {
    char* cur_line = getline(f_sam);

    if(cur_line == NULL)
    {
      break;
    }

    // If this is a comment line, skip it.
    if(cur_line[0] == '@')
    {
      delete [] cur_line;
      continue;
    }

    if(sscanf(cur_line, "%[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %[^\t] %*[^\t] %*[^\t] %*[^\t] %[^\t] %[^\t]", read_id, flag_str, chrom, _chr_index_str, mapp_quality_str ,mapping_map_str, cur_fragment, phred_quality_str) == 8)
    {
      if(phred_score_base == -123)
      {
        int l_read = strlen(phred_quality_str);
        for(int i = 0; i < l_read; i++)
        {
          if(phred_quality_str[i] > 'J')
          {
            fprintf(stderr, "Phred+64 encoding @ %d. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
              cur_fragment, phred_quality_str);
            phred_score_base = ';';
            break;
          }
          else if(phred_quality_str[i] < ';')
          {
            fprintf(stderr, "Phred+33 encoding @ %d. read (%c):\n%s\n%s\n", n_processed_reads, phred_quality_str[i],
              cur_fragment, phred_quality_str);
            phred_score_base = '!';
            break;
          }
        } // i loop.
      }

      //fprintf(stderr, "Processing read @ %d parsed with %s:\n%s\n", chr_index, mapping_map_str, cur_fragment);
      // If the quality is not adequate, do not use this read.
      if(atoi(mapp_quality_str) < min_mapp_qual)
      {
        n_low_quality_reads++;
        delete [] cur_line;
        continue;
      }

      // Make sure the normalized chromosome ids match.
      normalize_chr_id(chrom);

      int flag = atoi(flag_str);

      int i_chr = t_string::get_i_str(chr_ids, chrom);

      // If we do not have the chromosome, do not process.
      if(i_chr == (int)chr_ids->size())
      {
        n_unmapped_reads++;
        delete [] cur_line;
        continue;
      }

      int _chr_index = atoi(_chr_index_str);

      // Translate the 0 based index in SAM file to codebase's indexing, which is 1 based inclusive.
      chr_index = translate_coord(_chr_index, SAM_COORDS::start_base, CODEBASE_COORDS::start_base);

      // Sanity check. Is this fragment mapped?
      if(flag & 0x04)
      {
        // The read is not mapping.
      }
      else
      {
        n_processed_reads++;
          
        if(n_processed_reads % 1000000 == 0)
        {
          fprintf(stderr, "Processing %d. read             \r", n_processed_reads);
        }

if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
        fprintf(stderr, "Processing:\n%s\n%s\n", cur_fragment, phred_quality_str);
}

        int i_mapp_map = 0;
        //t_string* cur_entry_length_str = new t_string();
        bool is_matching = false;
        char entry_type_char;

        // Parse the cigar string to get the fragments.
        bool is_read_spliced = false;
        bool mapping_map_str_valid = validate_mapping_map_str(mapping_map_str, is_read_spliced);

        int read_nuc_index = 0;

        // Check if the mapping map string has splicing information, if it does and there extension length is not 0, 
        while(mapping_map_str_valid && 
          mapping_map_str[i_mapp_map] != 0)
        {
          int l_cur_entry;
          get_next_entry_per_mapp_map_string(mapping_map_str,
                            i_mapp_map, 
                            is_matching,
                            l_cur_entry,
                            entry_type_char);

          // If this block is matching, update the pileup.
          if(is_matching)
          {
            // Update the counts for the nucleotides.
            int cur_read_nuc_i = read_nuc_index;
            for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
            {
              int cur_nuc_num = per_nuc_val[(int)cur_fragment[cur_read_nuc_i]];

              //if(numerized_sequence_signal[cur_genome_i] > 0 &&
              //  nuc_2_num(cur_fragment[cur_read_nuc_i]) != (numerized_sequence_signal[cur_genome_i]))
              //{
              //  fprintf(stderr, "Found %c->%c @ %d\n", num_2_nuc(numerized_sequence_signal[cur_genome_i]), cur_fragment[cur_read_nuc_i], cur_genome_i);
              //}

              if(cur_nuc_num < 4)
              {
                // Update the count: The allelic counts must be checked for bounds.
                if(per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i] < max_n_alleles_per_posn-5)
                {
                  if(phred_quality_str[cur_read_nuc_i] - phred_score_base > min_phred_qual)
                  {
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
                    fprintf(stderr, "Adding %d: %c, %c (%d)\n", 
                        cur_read_nuc_i, 
                        cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
                        (int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
}

                    // Phred quality holds, update the pileup position.
                    per_chrom_nuc_count_per_allele[i_chr][cur_nuc_num][cur_genome_i]++;
                  } // phred qual pass check.
                  else
                  {
if(__DUMP_MAPPED_READ_TOOLS_MSGS__)
{
                    fprintf(stderr, "Skipping %d: %c, %c (%d)\n", 
                      cur_read_nuc_i, 
                      cur_fragment[cur_read_nuc_i], phred_quality_str[cur_read_nuc_i], 
                      (int)(phred_quality_str[cur_read_nuc_i]-phred_score_base));
}
                  } // phred qual nopass check.
                } // max_n_alleles_per_posn check.
              } // char < 4 check.

              // Update the read's nucleotide index.
              cur_read_nuc_i++;
            } // i_cur loop.
          } // genomic match check.
          else if(entry_type_char == 'D')
          {
            // Deletion from the reference: Update the 4th entry: Add all of these entries as deletions.
            for(int cur_genome_i = chr_index; cur_genome_i <= chr_index+l_cur_entry-1; cur_genome_i++)
            {
              per_chrom_nuc_count_per_allele[i_chr][4][cur_genome_i]++;
            } // cur_genome_i loop.
          }
          else if(entry_type_char == 'I')
          {
            // Insertion to the reference: This is included to one position.
            per_chrom_nuc_count_per_allele[i_chr][4][chr_index]++;
          }
          // Update the base for the current entry.
          if(check_genome_index_update_per_CIGAR_entry(entry_type_char))
          {
            chr_index += l_cur_entry;
          }

          // Update the base for the current read if requested.
          if(check_read_nuc_index_update_per_CIGAR_entry(entry_type_char))
          {
            read_nuc_index += l_cur_entry;
          }
        } // mapping map string processing loop.

        //delete(cur_entry_length_str);
      } // mapping check for the current read.
    } // SAM read line parse check.
    else
    {
      fprintf(stderr, "Could not parse %s\n", cur_line);
      exit(0);
    }

    delete [] cur_line;
  } // file reading loop.

  fprintf(stderr, "Finished reading.\n");

  fclose(f_sam);

  // Dump the counts per position.
  for(int i_chr = 0; i_chr < (int)chr_ids->size(); i_chr++)
  {
    char cur_allele_fp[1000];
    sprintf(cur_allele_fp, "%s/%s_allele_counts.bin", op_dir, chr_ids->at(i_chr));

    // compress and dump the current file.
    compress_nucleotide_pileup_track(per_chrom_nuc_count_per_allele[i_chr], chr_lengths->at(i_chr), cur_allele_fp);
  } // i_chr loop.

  //char analysis_summary_str[1000];
  //sprintf(analysis_summary_str, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality, n_unmapped_reads);
  //t_config_params::copy_analysis_summary_string(analysis_summary_str);

  //// Write the summary file.
  //f_summary = open_f(summary_fp, "w");
  //fprintf(f_summary, "Processed %d reads, %d low quality reads, %d unmapped reads.\n", n_processed_reads, n_low_quality_reads, n_unmapped_reads);
  //fclose(f_summary);
}
char* load_binary_sequence_file(char* bin_seq_fp, int& l_seq)
{
	FILE* f_bin_seq = open_f(bin_seq_fp, "rb");
	int l_bin_seq = 0;

	// Read the length.
	fread(&l_bin_seq, sizeof(int), 1, f_bin_seq);

	l_seq = l_bin_seq;

	// Read the sequence.
	char* bin_seq = new char[l_bin_seq + 2];
	fread(bin_seq, sizeof(char), l_bin_seq + 1, f_bin_seq);
	fclose(f_bin_seq);

	return(bin_seq);
}
void delete_pileup(unsigned short** loaded_pileup)
{
	//unsigned short** loaded_pileup = new unsigned short*[5];
	for (int i_allele = 0; i_allele < 5; i_allele++)
	{
		delete[] loaded_pileup[i_allele];
	} // i_allele loop.

	delete[] loaded_pileup;
}
int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "USAGE: %s [Options] [Arguments]\n\
Options:\n\
	-generate_compressed_pileup_per_SAM [SAM file path (\"stdin\" for piping from samtools)] [Chromosome Ids/lengths file path] [Output directory] [Minimum mapping quality] [Minimum base quality]\n\
	-get_SNVs_per_pileup [chromosome info file path] [Pileup directory] [Binary sequences directory] [Min coverage per SNV (20)] [Min MAF covg per SNV (4)] [Min MAF (0.2)] [Output file path]\n", argv[0]);
		exit(0);
	}

	if (t_string::compare_strings(argv[1], "-generate_compressed_pileup_per_SAM"))
	{
		if (argc != 7)
		{
			fprintf(stderr, "USAGE: %s -generate_compressed_pileup_per_SAM [SAM file path (\"stdin\" for piping from samtools)] [Chromosome Ids/lengths file path] [Output directory] [Minimum mapping quality] [Minimum base quality]\n", argv[0]);
			exit(0);
		}

		char* sam_fp = argv[2];
		char* chr_ids_lengths_list_fp = argv[3];
		char* op_dir = argv[4];
		int min_mapp_qual = atoi(argv[5]);
		int min_base_qual = atoi(argv[6]);

		vector<char*>* chr_ids = new vector<char*>();
		vector<int>* chr_lengths = new vector<int>();
		load_chromosome_lengths_per_tabbed_file(chr_ids_lengths_list_fp, chr_ids, chr_lengths);

		int n_processed_reads = 0;
		dump_nucleotide_pileup_per_SAM_file(sam_fp, chr_ids, chr_lengths, op_dir, min_mapp_qual, min_base_qual, n_processed_reads);
	} // -generate_compressed_pileup_per_SAM	
	else if (t_string::compare_strings(argv[1], "-get_SNVs_per_pileup"))
	{
		if (argc != 9)
		{
			fprintf(stderr, "USAGE: %s -get_SNVs_per_pileup [chromosome info file path] [Pileup directory] [Binary sequences directory] [Min coverage per SNV (20)] [Min MAF covg per SNV (4)] [Min MAF (0.2)] [Output file path]\n", argv[0]);
			exit(0);
		}

		char* chr_info_fp = argv[2];
		char* pileup_dir = argv[3];
		char* bin_seq_dir = argv[4];
		int min_covg = atoi(argv[5]);
		int min_alternate_covg = atoi(argv[6]);
		double min_alternate_freq = atof(argv[7]);
		char* op_fp = argv[8];

		fprintf(stderr, "Calling SNVs with minimum %d total and %d alternate read coverage and %.3f min alternate frequency.\n", min_covg, min_alternate_covg, min_alternate_freq);

		vector<char*>* chr_ids = new vector<char*>();
		vector<int>* chr_lengths = new vector<int>();

		int* nuc_2_num_array = get_per_char_as_nuc_2_num_coding_array();
		load_chromosome_lengths_per_tabbed_file(chr_info_fp, chr_ids, chr_lengths);

		FILE* f_op = open_f(op_fp, "w");
		for (int i_chr = 0; i_chr < chr_ids->size(); i_chr++)
		{
			fprintf(stderr, "Procesing %s\n", chr_ids->at(i_chr));

			int l_cur_chr_seq = 0;
			char cur_chr_bin_fp[1000];
			sprintf(cur_chr_bin_fp, "%s/%s.bin", bin_seq_dir, chr_ids->at(i_chr));
			char* cur_chr_seq = load_binary_sequence_file(cur_chr_bin_fp, l_cur_chr_seq);

			char cur_chr_pileup_fp[1000];
			sprintf(cur_chr_pileup_fp, "%s/%s_allele_counts.bin", pileup_dir, chr_ids->at(i_chr));
			if (!check_file(cur_chr_pileup_fp))
			{
				sprintf(cur_chr_pileup_fp, "%s/%s_allele_counts.bin.gz", pileup_dir, chr_ids->at(i_chr));
				if (!check_file(cur_chr_pileup_fp))
				{
					fprintf(stderr, "Could not find %s.\n", cur_chr_pileup_fp);
					exit(0);
				}
			}

			int l_pileup = 0;
			unsigned short** cur_chr_pileup = load_compressed_pileups(cur_chr_pileup_fp, l_pileup);
			fprintf(stderr, "Loaded %d long pileup.\n", l_pileup);

			char per_num_nucs[] = "ACGTN";

			for (int i = 1; i <= l_pileup; i++)
			{
				int cur_total_covg = 0;
				int max_alternate_cnt = 0;
				int max_alt_allele_i = -1;
				for (int i_all = 0; i_all < 4; i_all++)
				{
					cur_total_covg += cur_chr_pileup[i_all][i];

					bool is_alternate_allele = (nuc_2_num_array[cur_chr_seq[i]] != i_all);
					if (is_alternate_allele &&
						cur_chr_pileup[i_all][i] > min_alternate_covg &&
						max_alternate_cnt < cur_chr_pileup[i_all][i])
					{
						max_alternate_cnt = cur_chr_pileup[i_all][i];
						max_alt_allele_i = i_all;
					}
				} // i_all loop.

				if (cur_total_covg < min_covg)
				{
					continue;
				}

				double cur_alternate_freq = (double)(max_alternate_cnt) / (double)(cur_total_covg);
				if (cur_alternate_freq < min_alternate_freq)
				{
					continue;
				}

				fprintf(f_op, "%s\t%d\t%c\t%c\t%d\t%d\n", chr_ids->at(i_chr), i,
					cur_chr_seq[i], per_num_nucs[max_alt_allele_i],
					max_alternate_cnt, cur_total_covg);
			} // i loop.

			  // Free memory.
			delete_pileup(cur_chr_pileup);

			// Free sequence memory.
			delete[] cur_chr_seq;
		} // i_chr loop.

		fclose(f_op);
	} // -get_SNVs_per_pileup option.
}


