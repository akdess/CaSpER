#ifndef __MAPPED_READ_FILE_INTERFACE__
#define __MAPPED_READ_FILE_INTERFACE__

#include <vector>
using namespace std;

const int MEG_BASE = 1000 * 1000;
const int K_BASE = 1000;

class t_string;
struct t_annot_region;

bool check_genome_index_update_per_CIGAR_entry(char entry_char);
bool check_read_nuc_index_update_per_CIGAR_entry(char entry_char);

void dump_nucleotide_pileup_per_SAM_file(char* sam_fp, vector<char*>* chr_ids, vector<int>* chr_lengths, char* op_dir, int min_mapp_qual, int min_phred_qual, int& n_processed_reads);
void compress_nucleotide_pileup_track(unsigned short** pileup, int l_sig, char* op_fp);
unsigned short** allocate_pileup(int l_sig);
void delete_pileup(unsigned short** loaded_pileup);
unsigned short** load_compressed_pileups(char* cur_comp_allele_fp, int& l_pileup);

bool validate_mapping_map_str(char* mapping_map_str, bool& is_read_spliced);
void get_next_entry_per_mapp_map_string(char* mapping_map_str,
										int& i_mapp_map, 
										bool& is_matching,
										int& l_cur_entry,
										char& entry_type_char);

#endif // __MAPPED_READ_FILE_INTERFACE__

