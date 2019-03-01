#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h> // For va_... functions.
#include "ansi_string.h"
#include <ctype.h>
#include <algorithm>

// Simple Thread-safe string library.
// DO NOT INCLUDE string.h!

t_string::t_string(char* string)
{
	// Allocate the memory.
	//this->obj_string = (char*)malloc(string_length(string) + 3);
	//this->obj_str_mem_size = string_length(string) + 1;

	this->string_buffer = new t_string_buffer();
	this->string_buffer->i_buff = 0;
	this->string_buffer->l_buffer_mem = IBS;
	this->string_buffer->string_buff = new char[IBS+3];

	this->copy(string);
}

// This is the memset from string library.
void t_string::set_byte_buffer(void* buffer, long int n_bytes_2_set, char val_2_set)
{
	uint8_t* byte_buff = (uint8_t*)buffer;
	for(long int i_8b = 0; i_8b < n_bytes_2_set; i_8b++)
	{
		byte_buff[i_8b] = val_2_set;
	} // i_char loop.

	//uint64_t* large_buff = (uint64_t*)buffer;
	//long int n_64b_vals = n_bytes_2_set / 8;
	//int n_8b_vals = n_bytes_2_set % 8;
	//fprintf(stderr, "%d 64b vals, %d 8b vals.\n", n_64b_vals, n_8b_vals);
	//uint64_t large_val_2_set = (val_2_set << 56) + 
	//							(val_2_set << 48) + 
	//							(val_2_set << 40) + 
	//							(val_2_set << 32) + 
	//							(val_2_set << 24) + 
	//							(val_2_set << 16) + 
	//							(val_2_set << 8) + 
	//							val_2_set;

	//for(long int i_64b = 0; i_64b < n_64b_vals; i_64b++)
	//{
	//	large_buff[i_64b] = large_val_2_set;
	//} // i_char loop.

	//uint8_t* byte_buff = (uint8_t*)buffer;
	//for(int i_8b = n_64b_vals*8; i_8b < n_64b_vals*8+n_8b_vals; i_8b++)
	//{
	//	byte_buff[i_8b] = val_2_set;
	//} // i_char loop.
}

t_string::t_string(t_string* string)
{
	// Allocate the memory.
	//this->obj_string = (char*)malloc(string_length(string) + 3);
	//this->obj_str_mem_size = string_length(string) + 1;
	this->string_buffer = new t_string_buffer();
	this->string_buffer->i_buff = 0;
	this->string_buffer->l_buffer_mem = IBS;
	this->string_buffer->string_buff = new char[IBS+3];

	this->copy(string);
}

t_string::t_string()
{
	//this->obj_string = (char*)malloc(IBS + 3);
	//this->obj_str_mem_size = IBS;
	this->string_buffer = new t_string_buffer();
	this->string_buffer->i_buff = 0;
	this->string_buffer->l_buffer_mem = IBS;
	this->string_buffer->string_buff = new char[IBS+3];

	// Initialize the string.
	this->empty();
}

bool t_string::sort_strings_per_prefix(char* str1, char* str2)
{
	int l1 = t_string::string_length(str1);
	int l2 = t_string::string_length(str2);
	int max_l12 = (l1 > l2)?(l1):(l2);
	for(int i = 0; i < max_l12; i++)
	{
		// First string is shorter and the prefixes match exactly, put str1 frist.
		if(str1[i] == 0)
		{
			return(true);
		}

		// The second string is shorter, it should come before first.
		if(str2[i] == 0)
		{
			return(false);
		}

		if(str1[i] < str2[i])
		{
			return(true);
		}
		else if(str1[i] > str2[i])
		{
			return(false);
		}
	} // i loop.

	// Strings are exactly the same, return false to obey strict weak ordering.
	return(false);
}

bool t_string::sort_strings(char* str1, char* str2)
{
	int l1 = t_string::string_length(str1);
	int l2 = t_string::string_length(str2);
	if(l1 < l2)
	{
		return(true);
	}
	else if(l1 > l2)
	{
		return(false);
	}
	else
	{
		for(int i = 0; i < l1; i++)
		{
			if(str1[i] < str2[i])
			{
				return(true);
			}
			else if(str1[i] > str2[i])
			{
				return(false);
			}
		} // i loop.

		// Strings are exactly the same, return false to obey strict weak ordering.
		return(false);
	}
}

int t_string::get_matching_prefix_length(char* str1, char* str2)
{
	int l1 = string_length(str1);
	int l2 = string_length(str2);
	int l_str_2_match = (l1 < l2)?(l1):(l2);
	for(int i = 0; i < l_str_2_match; i++)
	{
		if(str1[i] != str2[i])
		{
			return(i);
		}
	} // i loop.

	return(l_str_2_match);
}

bool t_string::compare_strings_per_total_prefix(char* str1, char* str2)
{
	int l1 = t_string::string_length(str1);
	int l2 = t_string::string_length(str2);
	int l_match_prefix = t_string::get_matching_prefix_length(str1, str2);

	if(l1 == l_match_prefix ||
		l2 == l_match_prefix)
	{
		return(true);
	}

	return(false);
}


vector<int>* t_string::get_string_sorting_idx(vector<char*>* string_list)
{
	vector<t_str_node*>* str_nodes = new vector<t_str_node*>();
	for(int i_s = 0; i_s < (int)string_list->size(); i_s++)
	{
		t_str_node* new_node = new t_str_node();
		new_node->i = i_s;
		new_node->str = string_list->at(i_s);
		str_nodes->push_back(new_node);
	} // i_s loop.

	// Sort the nodes.
	sort(str_nodes->begin(), str_nodes->end(), t_string::sort_str_nodes);	

	// Fill the sorting idx.
	vector<int>* string_sort_idx = new vector<int>();
	for(int i_s = 0; i_s < (int)str_nodes->size(); i_s++)
	{
		string_sort_idx->push_back(str_nodes->at(i_s)->i);
		delete(str_nodes->at(i_s));
	} // i_s loop.

	delete(str_nodes);

	return(string_sort_idx);
}

bool t_string::sort_str_nodes(t_str_node* node1, t_str_node* node2)
{
	return(t_string::sort_strings(node1->str, node2->str));
}

/*
Match strings, otherwise return the first string in the list that is just before the query.
If the query is before first string in the list, return 0th string. 
*/
int t_string::fast_search_string_per_prefix(char* query, vector<char*>* string_list, int i, int j)
{
	// Sort the string list.
	int i_mid = (i + j) / 2;

	if(i_mid == i ||
		i_mid == j)
	{
		return(i);
	}
	
	if(t_string::compare_strings(string_list->at(i_mid), query))
	{
		while(i_mid < (int)string_list->size() &&
			t_string::sort_strings_per_prefix(string_list->at(i_mid), query))
		{
			i_mid++;
		}

		// If the above loop hit the end, move back once.
		if(i_mid == (int)string_list->size())
		{
			i_mid--;
		}

		// If the above loop moved beyond the query move back once, which is the exact element just before query.
		if(i_mid > 0 &&
			!t_string::sort_strings_per_prefix(string_list->at(i_mid), query))
		{
			i_mid--;
		}

		return(i_mid);
	}
	else if(t_string::sort_strings_per_prefix(query, string_list->at(i_mid)))
	{
		return(t_string::fast_search_string_per_prefix(query, string_list, i, i_mid));
	}
	else if(!t_string::sort_strings_per_prefix(query, string_list->at(i_mid)))
	{
		return(t_string::fast_search_string_per_prefix(query, string_list, i_mid, j));
	}
	else
	{
		fprintf(stderr, "Fast searching failed @ %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
}

int t_string::fast_search_string(char* query, vector<char*>* string_list, int i, int j)
{
	// Sort the string list.
	int i_mid = (i + j) / 2;

	if(i_mid == i ||
		i_mid == j)
	{
		return(i);
	}
	
	if(t_string::compare_strings(string_list->at(i_mid), query))
	{
		while(i_mid < (int)string_list->size() &&
			t_string::sort_strings(string_list->at(i_mid), query))
		{
			i_mid++;
		}

		// If the above loop hit the end, move back once.
		if(i_mid == (int)string_list->size())
		{
			i_mid--;
		}

		// If the above loop moved beyond the query move back once, which is the exact element just before query.
		if(i_mid > 0 &&
			!sort_strings(string_list->at(i_mid), query))
		{
			i_mid--;
		}

		return(i_mid);
	}
	else if(sort_strings(query, string_list->at(i_mid)))
	{
		return(t_string::fast_search_string(query, string_list, i, i_mid));
	}
	else if(!sort_strings(query, string_list->at(i_mid)))
	{
		return(t_string::fast_search_string(query, string_list, i_mid, j));
	}
	else
	{
		fprintf(stderr, "Fast searching failed @ %s(%d)\n", __FILE__, __LINE__);
		exit(0);
	}
}

int t_string::fast_count_non_empty_tokens(char* str, char* delims)
{
	int n_delims = t_string::string_length(delims);
	int n_toks = 0;
	int l_cur_tok = 0;
	int i = 0;
	while(str[i] != 0)
	{
		//bool found_delim = false;
		for(int i_d = 0; i_d < n_delims; i_d++)
		{
			if(str[i] == delims[i_d])
			{
				if(l_cur_tok > 0)
				{
					n_toks++;
					l_cur_tok = 0;
				}
				else
				{
					// This is an empty token, do not count it.
				}
			}
			else
			{
				// This is a non-delim char, update the current token length.
				l_cur_tok++;
				break;
			}
		} // i_d loop.

		i++;
	} // string loop.

	if(l_cur_tok > 0)
	{
		n_toks++;
	}

	return(n_toks);
}

char* t_string::copy_me_str(const char* str)
{
	if(str == NULL)
	{
		return(NULL);
	}

	char* new_str = new char[t_string::string_length(str) + 2];
	t_string::copy(new_str, str);
	return(new_str);
}

// Following function is useful for automatically adding a string to a list which does not include it. The returned index is the index of the string in the list.
int t_string::get_add_i_str(vector<char*>* strs, char* str)
{
	int i_chr = t_string::get_i_str(strs, str);

	// Add to the list.
	if(i_chr == (int)strs->size())
	{
		char* new_chr_id = t_string::copy_me_str(str);
		strs->push_back(new_chr_id);

		i_chr = t_string::get_i_str(strs, str);

		if(i_chr == (int)strs->size())
		{
			fprintf(stderr, "Failed to insert %s to the list.\n", str);
			return((int)strs->size());
		}
	}

	return(i_chr);
}

int t_string::get_i_str(vector<t_string*>* strs, char* str)
{
	for(int i_chr = 0; i_chr < (int)strs->size(); i_chr++)
	{
		if(t_string::compare_strings(strs->at(i_chr)->str(), str))
		{
			return(i_chr);
		}
	} // i_chr loop.

	//printf("Could not find chromosome id %s @ %s(%d), adding it to the list.\n", chrom_id, __FILE__, __LINE__);
	return((int)strs->size());
}

int t_string::get_i_str_ci(vector<char*>* strs, char* str)
{
	for(int i_chr = 0; i_chr < (int)strs->size(); i_chr++)
	{
		if(t_string::compare_strings_ci(strs->at(i_chr), str))
		{
			return(i_chr);
		}
	} // i_chr loop.

	//printf("Could not find chromosome id %s @ %s(%d), adding it to the list.\n", chrom_id, __FILE__, __LINE__);
	return((int)strs->size());
}

int t_string::get_i_str(vector<char*>* strs, char* str)
{
	for(int i_chr = 0; i_chr < (int)strs->size(); i_chr++)
	{
		if(t_string::compare_strings(strs->at(i_chr), str))
		{
			return(i_chr);
		}
	} // i_chr loop.

	//printf("Could not find chromosome id %s @ %s(%d), adding it to the list.\n", chrom_id, __FILE__, __LINE__);
	return((int)strs->size());
	////exit(0);
	//getc(stdin);

	//char* new_id = new char[strlen(chrom_id) + 2];
	//strcpy(new_id, chrom_id);
	//chr_ids->push_back(new_id);

	//return(t_string::get_i_str(chr_ids, chrom_id));
}

t_string::~t_string()
{
	//printf("Freeing a string object.\n");
	//free(this->obj_string);
	delete [] this->string_buffer->string_buff;
	delete this->string_buffer;
}

vector<int>* t_string::get_integers_in_string(char* str)
{
	vector<int>* int_vector = new vector<int>();
	int l_str = t_string::string_length(str);

	// The number string.
	t_string* cur_int_str = new t_string();

	for(int char_i = 0; char_i < l_str; char_i++)
	{
		if(str[char_i] >= '0' && str[char_i] <= '9')
		{
			// Add this character.
			cur_int_str->concat_char(str[char_i]);
		}
		else
		{
			if(cur_int_str->length() > 0)
			{
				int new_num = str2num(cur_int_str, 10);

				// Add the new number
				int_vector->push_back(new_num); 

				cur_int_str->empty();
			}
		}
	} // char index.

	// Process the last string if it is a number.
	if(cur_int_str->length() > 0)
	{
		int new_num = str2num(cur_int_str, 10);

		// Add the new number
		int_vector->push_back(new_num); 

		cur_int_str->empty();
	}

	// Free integer storing string.
	delete(cur_int_str);

	return(int_vector);
}

vector<int>* t_string::get_integers_in_string()
{
	vector<int>* int_vector = new vector<int>();

	// The number string.
	t_string* cur_int_str = new t_string();

	int l_str = this->length();
	for(int char_i = 0; char_i < l_str; char_i++)
	{
		if(this->x(char_i) >= '0' && this->x(char_i) <= '9')
		{
			// Add this character.
			cur_int_str->concat_char(this->x(char_i));
		}
		else
		{
			if(cur_int_str->length() > 0)
			{
				int new_num = str2num(cur_int_str, 10);

				// Add the new number
				int_vector->push_back(new_num); 

				cur_int_str->empty();
			}
		}
	} // char index.

	// Process the last string if it is a number.
	if(cur_int_str->length() > 0)
	{
		int new_num = str2num(cur_int_str, 10);

		// Add the new number
		int_vector->push_back(new_num); 

		cur_int_str->empty();
	}

	// Free integer storing string.
	delete(cur_int_str);

	return(int_vector);
}

void t_string::copy(char* dest_string, const char* src_string)
{
	unsigned int l_src_string = string_length(src_string);
	
	// Do the actual the copying, this loop copies the null character!
	for(unsigned int i = 0; i <= l_src_string; i++)
	{
		dest_string[i] = src_string[i];
	}
	dest_string[l_src_string] = 0;
}

/*
Copy does not require an initialized string in the object's buffer.
*/
void t_string::copy(const char* string)
{
	//fprintf(stderr, "Copying %s\n", string);

	//unsigned int l_src_string = string_length(string);
	////if(this->obj_str_mem_size < l_src_string + 3)
	//if(this->string_buffer->l_buffer_mem < l_src_string + 3)
	//{
	//	//free(this->obj_string);
	//	delete [] this->string_buffer->string_buff;
	//	this->string_buffer->string_buff = new char [l_src_string + 3];
	//	this->string_buffer->i_buff = l_src_string;
	//	this->string_buffer->l_buffer_mem = l_src_string + 3;

	//	//this->obj_str_mem_size = l_src_string + 3;
	//	//this->obj_string = (char*)malloc(l_src_string + 3);
	//	//this->obj_str_mem_size = l_src_string + 3;
	//}
	//
	//// Do the actual the copying, this loop copies the null character!
	//for(unsigned int i = 0; i <= l_src_string; i++)
	int i = 0;

	while(string[i] != 0)
	{
		//this->x(i) = string[i];
		this->concat_char(string[i]);
		i++;
	}

	//// Set the end of the string.
	//this->x(i) = 0;
}

char* t_string::get_alphanumeric(char* string)
{
	int l_str = t_string::string_length(string);
	char* alphanum_string = new char[l_str + 2];
	//memset(alphanum_string, 0, l_str);
	t_string::set_byte_buffer(alphanum_string, l_str, 0);
	int i_an = 0;
	for(int i = 0; i < l_str; i++)
	{
		// Update if this is an alphanumeric character.
		if((string[i] >= 'a' && string[i] <= 'z') ||
			(string[i] >= 'A' && string[i] <= 'Z') ||
			(string[i] >= '0' && string[i] <= '9'))
		{
			alphanum_string[i_an] = string[i];
			i_an++;
		}
	} // i loop.

	return(alphanum_string);
}

void t_string::copy(t_string* string)
{
	this->copy(string->string_buffer->string_buff);
}

// Static string library functions: These are thread safe functions.
int t_string::string_length(const char* string)
{
	unsigned int str_length = 0;
	while(string[str_length] != 0)
	{
		str_length++;
	}

	return(str_length);
}

// Static string library functions: These are thread safe functions.
int t_string::string_length(t_string* string)
{
	//return(string->length());
	return(string->string_buffer->i_buff);
}

int t_string::length()
{
	return(string_length(this));
}

// empty the string.
void t_string::empty()
{
	this->string_buffer->i_buff = 0;
	this->x(0) = 0;
}

// Getter function.
char& t_string::x(int i)
{
	//return(this->obj_string[i]);
	return(this->string_buffer->string_buff[i]);
}

char* t_string::str()
{
	//return(this->obj_string);
	return(this->string_buffer->string_buff);
}

char* t_string::char_2_str(const char char_val)
{
	char* char_str = new char[2];
	char_str[0] = char_val;
	char_str[1] = 0;

	return char_str;
}

char* t_string::substring(const char* str, int i, int j)
{
	//t_string* sub_str = new t_string();
	char* sub_str = NULL;

	int l_str = t_string::string_length(str);

	// Validity checks on the arguments.
	if(i > j ||
		i > l_str)
	{
		return(NULL);
	}
	else if(j >= l_str)
	{
		j = l_str-1;
	}

	sub_str = new char[j-i+1+2];
	sub_str[0] = 0; // Make this an empty string.

	// Copy the substring.
	int ip = i;
	int i_sub_str = 0;
	while(ip <= j)
	{
		sub_str[i_sub_str] = str[ip];

		i_sub_str++;
		ip++;
	}

	// End the string.
	sub_str[i_sub_str] = 0;

	return(sub_str);
}

// Return the substring at substring starting at i and ending at j, inclusive.
char* t_string::substring(int i, int j)
{
	//t_string* sub_str = new t_string();
	char* sub_str;

	// Validity checks on the arguments.
	if(i > j ||
		i > this->length() || 
		j > this->length())
	{
		return(NULL);
	}
	else
	{
		sub_str = (char*)malloc(sizeof(char) * (j-i+1+2));
		sub_str[0] = 0; // Make this an empty string.
	}

	// Copy the substring.
	int ip = i;
	int i_sub_str = 0;
	while(ip <= j)
	{
		sub_str[i_sub_str] = this->x(ip);

		i_sub_str++;
		ip++;
	}

	// End the string.
	sub_str[i_sub_str] = 0;

	return(sub_str);
}

// This is a simple fuzzy string comparison.
char* t_string::get_longest_matching_substring(char* string1, char* string2)
{
	return(NULL);
}

bool t_string::is_empty(char* string)
{
	int l_str = t_string::string_length(string);

	for(int i = 0; i < l_str; i++)
	{
		// 
		if(string[i] != ' ' && 
			string[i] != '\t' &&
			string[i] != '\n')
		{
			return(false);
		}
	}

	return(true);
}

enum{STR_STATE_PREFIX, STR_STATE_WHOLE, STR_STATE_DECIMAL, STR_STATE_EXPO};
bool t_string::is_number(char* string)
{
	int cur_state = STR_STATE_PREFIX;

	for(int i = 0; ; i++)
	{
		if(string[i] == 0)
		{
			break;
		}

		if(string[i] == '+' || 
			string[i] == '-')
		{
			if(cur_state == STR_STATE_EXPO)
			{
			}
			else
			{
				cur_state = STR_STATE_WHOLE;
			}
		}
		else if(string[i] == '.')
		{
			// We can come to this state from whole or prefix state, not from decimal state.
			if(cur_state == STR_STATE_DECIMAL)
			{
				return(false);
			}

			cur_state = STR_STATE_DECIMAL;
		}
		else if(string[i] == ',')
		{
			// This is not allowed in decimal part of a number.
			if(cur_state == STR_STATE_DECIMAL ||
				cur_state == STR_STATE_PREFIX)
			{
				return(false);
			}
		}
		else if(string[i] >= '0' && string[i] <= '9')
		{	
			// Nothing changes.
		}
		else if(string[i] == 'e')
		{
			// Cannot get into this state again.
			if(cur_state == STR_STATE_EXPO)
			{
				return(false);
			}

			// If the state was prefix, move to whole.
			cur_state = STR_STATE_EXPO;

			if(string[i+1] == 0)
			{
				break;
			}
			
			// Otherwise nothing changes.
		}
		else
		{
			// Any other character is not acceptable as a number.
			return(false);
		}
	} // i loop.

	return(true);
}

/*
bool t_string::_is_number(char* string)
{
	int l_str = t_string::string_length(string);

	bool found_whole_part = false;
	bool found_decimal_part = false;

	for(int i = 0; i < l_str; i++)
	{
		if(string[i] == '.')
		{
			// The decimal place was already found, return false.
			if(found_decimal_part)
			{
				return(false);
			}

			// Found the decimal place.
			found_decimal_part = true;
		}
		else if(string[i] == ',')
		{
		}
		else if(string[i] >= '0' && string[i] <= '9')
		{
			// Found the whole place.
			if(!found_whole_part)
			{
				found_whole_part = true;
			}
		}
		else
		{
			// Any other character is not acceptable as a number.
			return(false);
		}
	} // i loop.

	return(true);
}
*/
/*
Replace the avoided chars in the list with a char. Note that setting char_to_replace to 0 removes the characters in the avoid list 
from the original string.
*/
void t_string::replace_avoid_list(char* str, char start_char, char end_char, char char_to_replace)
{
	t_string* avoided_char_list_str = new t_string();
	for(char cur_char = start_char; cur_char <= end_char; cur_char++)
	{
		avoided_char_list_str->concat_char(cur_char);
	} // cur_char loop.

	t_string::replace_avoid_list(str, avoided_char_list_str->str(), char_to_replace);
}

void t_string::replace_avoid_list(char* str, char* avoided_char_list, char char_to_replace)
{
	char* temp_str = new char[string_length(str) + 3];
	//memset(temp_str, 0, string_length(str) + 3);
	t_string::set_byte_buffer(temp_str, string_length(str) + 3, 0);
	int i_tmp_str = 0;
	for(int i_str = 0; i_str < string_length(str); i_str++)
	{
		// Look for avoided chars. if found, replace.
		bool found_avoided_char = false;
		for(int i_avoided = 0; i_avoided < string_length(avoided_char_list); i_avoided++)
		{
			if(str[i_str] == avoided_char_list[i_avoided])
			{
				found_avoided_char = true;
			}
		} // i_avoided loop

		if(found_avoided_char)
		{
			if(char_to_replace == 0)
			{
				// Do not put anything in the temporary string.
			}
			else
			{
				temp_str[i_tmp_str] = char_to_replace;
				i_tmp_str++;
			}
		}
		else
		{
			temp_str[i_tmp_str] = str[i_str];
			i_tmp_str++;
		}
	} // i_str loop

	copy(str, temp_str);
	delete [] temp_str;
}

vector<void*>* t_string::pool_obj_list(vector<vector<void*>*>* obj_lists_list)
{
	vector<void*>* pooled_list = new vector<void*>();

	for(int i_l = 0; i_l < (int)obj_lists_list->size(); i_l++)
	{
		pooled_list->insert(pooled_list->end(), obj_lists_list->at(i_l)->begin(), obj_lists_list->at(i_l)->end());
	} // i_l loop.

	return(pooled_list);
}

vector<void*>* t_string::interleave_obj_list(vector<vector<void*>*>* obj_lists_list)
{
	vector<void*>* interleaved_list = new vector<void*>();

	int* indices_per_list = new int[(int)obj_lists_list->size()];

	// At the begining all the indices are 0.
	for(int i_i = 0; i_i < (int)obj_lists_list->size(); i_i++)
	{

		indices_per_list[i_i] = 0;
	} // i_i loop.

	while(1)
	{
		bool added = false;

		for(int i_l = 0; i_l < (int)obj_lists_list->size(); i_l++)
		{
			// If the index for this list did not pass the end of the list, add the object to the list.
			if(indices_per_list[i_l] < (int)obj_lists_list->at(i_l)->size())
			{
				interleaved_list->push_back(obj_lists_list->at(i_l)->at(indices_per_list[i_l]));
				indices_per_list[i_l]++;
				added = true;
			}
		} // i_l loop.

		if(!added)
		{
			return(interleaved_list);
		}
	} // infinite loop exits when there are no more obj's to add.
}

vector<vector<void*>*>* t_string::split_object_list(vector<void*>* obj_list, int n_batches)
{
	vector<vector<void*>*>* batches = new vector<vector<void*>*>();

	int n_objs_per_batch = (int)obj_list->size() / n_batches;
	vector<void*>* current_batch = new vector<void*>();
	batches->push_back(current_batch);
	for(int i_str = 0; i_str < (int)obj_list->size(); i_str++)
	{
		if((int)current_batch->size() == n_objs_per_batch)
		{
			current_batch = new vector<void*>();
			batches->push_back(current_batch);
		}

		current_batch->push_back(obj_list->at(i_str));
	} // i_str loop.

	return(batches);
}

vector<vector<char*>*>* t_string::split_string_list(vector<char*>* string_list, int n_batches)
{
	vector<vector<char*>*>* batches = new vector<vector<char*>*>();

	int n_string_per_batch = (int)string_list->size() / n_batches;
	vector<char*>* current_batch = new vector<char*>();
	batches->push_back(current_batch);
	for(int i_str = 0; i_str < (int)string_list->size(); i_str++)
	{
		if((int)current_batch->size() == n_string_per_batch)
		{
			current_batch = new vector<char*>();
			batches->push_back(current_batch);
		}

		current_batch->push_back(string_list->at(i_str));
	} // i_str loop.

	return(batches);
}

vector<vector<int>*>* t_string::split_int_list(vector<int>* int_list, int n_batches)
{
	vector<vector<int>*>* batches = new vector<vector<int>*>();

	int n_ints_per_batch = (int)int_list->size() / n_batches;
	vector<int>* current_batch = new vector<int>();
	batches->push_back(current_batch);
	for(int i_str = 0; i_str < (int)int_list->size(); i_str++)
	{
		if((int)current_batch->size() == n_ints_per_batch)
		{
			current_batch = new vector<int>();
			batches->push_back(current_batch);
		}

		current_batch->push_back(int_list->at(i_str));
	} // i_str loop.

	return(batches);
}

vector<char*>* t_string::get_unique_entries(vector<char*>* str_list)
{
	if(str_list == NULL)
	{
		return(NULL);
	}

	vector<char*>* unique_entries = new vector<char*>();
	for(int i_str = 0; i_str < (int)str_list->size(); i_str++)
	{
		int i_ent = t_string::get_i_str(unique_entries, str_list->at(i_str));
		if(i_ent == (int)unique_entries->size())
		{
			unique_entries->push_back(t_string::copy_me_str(str_list->at(i_str)));
		}
	} // i_str loop.

	return(unique_entries);
}

bool t_string::get_next_token(char* string, char* buffer, int l_buffer, const char* delimiter_list, int& i_cur_char)
{
	int i_buff_char = 0;

	// Initialize the string.
	buffer[i_buff_char] = 0;
	while(string[i_cur_char] != 0)
	{
		// Search for delimiters.
		int i_delim = 0;
		while(delimiter_list[i_delim] != 0)
		{
			if(string[i_cur_char] == delimiter_list[i_delim])
			{
				// Move to the next character right after the delimiter.
				i_cur_char++;

				buffer[i_buff_char] = 0;

				if(i_buff_char > 0)
				{	
					return true;
				}
				else
				{
					return false;
				}
			}

			i_delim++;
		} // delimiter list search.

		// Copy the current character to the buffer and update the indices.
		if(i_buff_char < l_buffer)
		{
			buffer[i_buff_char] = string[i_cur_char];
		}

		i_cur_char++;
		i_buff_char++;
	} // string parsing loop.

	buffer[i_buff_char] = 0;

	if(i_buff_char > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

// Tokenize the first n tokens from a most probably very long string. The idea is to not process the whole string and to copy it.
vector<t_string*>* t_string::get_first_n_tokens(char* string, int n_tokens, const char* delimiter_list, int& i_next_char)
{
	vector<t_string*>* first_tokens = new vector<t_string*>();
	
	// This is the pointer to the next character in the string.
	int i_s = 0;

	// Process the next character.
	t_string* next_token = new t_string();
	while((int)first_tokens->size() != n_tokens &&
		string[i_s] != 0)
	{
		bool is_delim = false;
		int n_delimiters = t_string::string_length(delimiter_list);
		for(int i_d = 0; 
			!is_delim && i_d < n_delimiters; 
			i_d++)
		{
			if(delimiter_list[i_d] == string[i_s])
			{
				is_delim = true;
			}
		} // i_d loop.

		// This is not a delimiter.
		if(!is_delim)
		{
			next_token->concat_char(string[i_s]);
		} // delimiter check.
		else
		{
			// If this token contains characters, add it.
			if(next_token->length() > 0)
			{
				// Add the next token.
				first_tokens->push_back(next_token);

				// Reallocate the next token.
				next_token = new t_string();
			}
			else
			{
				delete(next_token);
				next_token = new t_string();
			}
		}

		//fprintf(stderr, "Procesing %d: %d. tokens: %s\n", i_s, first_tokens->size(), next_token->str());
		i_s++;
	} // token addition loop.

	// If the next token was generated for the last string in the string, add it to the list of tokens.
	// This is not added because the token adding is done per delimiter character.
	if(next_token->length() > 0)
	{
		first_tokens->push_back(next_token);
	}
	else
	{
		delete(next_token);
	}

	// Skip over all the delimiters.	
	bool is_delim = true;
	while(is_delim && string[i_s] != 0)
	{
		is_delim = false;

		int n_delimiters = t_string::string_length(delimiter_list);
		for(int i_d = 0; 
			!is_delim && i_d < n_delimiters; 
			i_d++)
		{
			if(delimiter_list[i_d] == string[i_s])
			{
				is_delim = true;
			}
		} // i_d loop.

		if(is_delim)
		{
			i_s++;
		}
	} // delimiter skipping loop.

	// Next character is the index of next character in the string, which is a delimiter.
	i_next_char = i_s;

	return(first_tokens);
}

// Tokenizer: Tokenizes with respect to the characters in the delimiter list, which is a null terminated list.
t_string_tokens* t_string::tokenize_by_chars(const char* delimiter_list)
{
	//printf("String: %sFIN\n", this->str());
	t_string_tokens* token_vector = new vector<t_string*>();

	// Current token is nothing at the beginning.
	t_string* current_token = new t_string();

	// Go over the obj_string.
	int str_len = this->length();
	int n_delimiters = string_length(delimiter_list);
	for(int i_str = 0; i_str < str_len; i_str++)
	{
		//printf("Processing: %c\n", this->x(i_str));
		// This is set to true if i_str'th character is a delimiter.
		bool is_delimiter = false;

		// Check if current character in the string is a delimiter character.
		for(int i_del = 0; i_del < n_delimiters; i_del++)
		{
			if(this->x(i_str) == delimiter_list[i_del])
			{
				is_delimiter = true;

				// Push current token to vector (if it is appropriate for pushing.) and start a new token.
				// Push the previous token. Cannot push the next token since there is no way of knowing that
				// it is a valid token beforehand.
				if(current_token != NULL && current_token->length() != 0 && !t_string::is_empty(current_token->str()))
				{
					// Push this token to token vector.
					//printf("pushing %s\n", current_token->obj_string);
					token_vector->push_back(current_token);

					// Do not free current token now, allocate a new one.
					current_token = new t_string();
				}
				else // The token is not a good token (i.e., an empty token), free its memory and 
				{
					current_token->empty();
				}

				break; // Break from the delimiter char loop.
			}
		} // delimiter char loop.

		// If this character is not a delimiter, concatenate the character to the current_token.
		if(!is_delimiter)
		{
			current_token->concat_char(this->x(i_str));
		}

	} // string char loop.

	// The last token is added to the tokens vector after whole loop is finished.
	if(current_token != NULL && current_token->length() != 0)
	{
		// Push this token to token vector.
		token_vector->push_back(current_token);
		//printf("pushing %s\n", current_token->obj_string);
	}
	else	
	{
		delete(current_token);
	}

	return(token_vector);
}

// Fast string tokenizer that does not call string length counting. This should work fairly fast on 
// long data.
t_string_tokens* t_string::tokenize_by_chars(char* string, const char* delimiter_list, bool return_empty_tokens)
{
	//printf("String: %sFIN\n", this->str());
	t_string_tokens* token_vector = new vector<t_string*>();

	// Current token is nothing at the beginning.
	t_string* current_token = new t_string();

	// Go over the obj_string.
	//int str_len = t_string::string_length(string);
	int n_delimiters = string_length(delimiter_list);

	int i_str = 0;
	while(string[i_str] != 0)
	{
		//printf("Processing: %c\n", this->x(i_str));
		// This is set to true if i_str'th character is a delimiter.
		bool is_delimiter = false;

		// Check if current character in the string is a delimiter character.
		for(int i_del = 0; i_del < n_delimiters; i_del++)
		{
			//if(this->x(i_str) == delimiter_list[i_del])
			if(string[i_str] == delimiter_list[i_del])
			{
				is_delimiter = true;

				// Push current token to vector (if it is appropriate for pushing.) and start a new token.
				// Push the previous token. Cannot push the next token since there is no way of knowing that
				// it is a valid token beforehand.
				if(current_token != NULL && 
					(current_token->length() != 0 || return_empty_tokens))
				{
					// Push this token to token vector.
					//printf("pushing %s\n", current_token->obj_string);
					token_vector->push_back(current_token);

					// Do not free current token now, allocate a new one.
					current_token = new t_string();
				}
				else // The token is not a good token (i.e., an empty token), free its memory and 
				{
					current_token->empty();
				}

				break; // Break from the delimiter char loop.
			}
		} // delimiter char loop.

		// If this character is not a delimiter, concatenate the character to the current_token.
		if(!is_delimiter)
		{
			//current_token->concat_char(this->x(i_str));
			current_token->concat_char(string[i_str]);
		}

		i_str++;
	} // string char loop.

	// The last token is added to the tokens vector after whole loop is finished.
	if(current_token != NULL && 
		(current_token->length() != 0 || return_empty_tokens))
	{
		// Push this token to token vector.
		token_vector->push_back(current_token);
		//printf("pushing %s\n", current_token->obj_string);
	}
	else	
	{
		delete(current_token);
	}

	return(token_vector);
}

vector<char*>* t_string::copy_tokens_2_strs(t_string_tokens* toks)
{
	vector<char*>* strs = new vector<char*>();
	for (int i_t = 0; i_t < toks->size(); i_t++)
	{
		strs->push_back(toks->at(i_t)->str());
	} // i_t loop.

	return(strs);
}

// Fast string tokenizer that does not call string length counting. This should work fairly fast on 
// long data.
t_string_tokens* t_string::tokenize_by_chars(char* string, const char* delimiter_list)
{
	//printf("String: %sFIN\n", this->str());
	t_string_tokens* token_vector = new vector<t_string*>();

	// Current token is nothing at the beginning.
	t_string* current_token = new t_string();

	// Go over the obj_string.
	//int str_len = t_string::string_length(string);
	int n_delimiters = string_length(delimiter_list);

	int i_str = 0;
	while(string[i_str] != 0)
	{
		//printf("Processing: %c\n", this->x(i_str));
		// This is set to true if i_str'th character is a delimiter.
		bool is_delimiter = false;

		// Check if current character in the string is a delimiter character.
		for(int i_del = 0; i_del < n_delimiters; i_del++)
		{
			//if(this->x(i_str) == delimiter_list[i_del])
			if(string[i_str] == delimiter_list[i_del])
			{
				is_delimiter = true;

				// Push current token to vector (if it is appropriate for pushing.) and start a new token.
				// Push the previous token. Cannot push the next token since there is no way of knowing that
				// it is a valid token beforehand.
				if(current_token != NULL && current_token->length() != 0)
				{
					// Push this token to token vector.
					//printf("pushing %s\n", current_token->obj_string);
					token_vector->push_back(current_token);

					// Do not free current token now, allocate a new one.
					current_token = new t_string();
				}
				else // The token is not a good token (i.e., an empty token), free its memory and 
				{
					current_token->empty();
				}

				break; // Break from the delimiter char loop.
			}
		} // delimiter char loop.

		// If this character is not a delimiter, concatenate the character to the current_token.
		if(!is_delimiter)
		{
			//current_token->concat_char(this->x(i_str));
			current_token->concat_char(string[i_str]);
		}

		i_str++;
	} // string char loop.

	// The last token is added to the tokens vector after whole loop is finished.
	if(current_token != NULL && current_token->length() != 0)
	{
		// Push this token to token vector.
		token_vector->push_back(current_token);
		//printf("pushing %s\n", current_token->obj_string);
	}
	else	
	{
		delete(current_token);
	}

	return(token_vector);
}

bool t_string::is_balanced(char* left_pars, char* right_pars)
{
	return(t_string::is_balanced(this->str(), left_pars, right_pars));
}

bool t_string::is_balanced(char* str, char* left_pars, char* right_pars)
{
        int n_pars = t_string::string_length(left_pars);

        // This is the list of left parenthesis characters in the string.
        char* left_par_chars = (char*)malloc(sizeof(char) * (t_string::string_length(str) + 1));

        // Count the left and right parenthesis characters, they should be equaal in number, watch the structure along the way.
        int l_str = t_string::string_length(str);
        int l_left_par_chars = 0;
        left_par_chars[l_left_par_chars] = 0; // Initialize.
        for(int i = 0; i < l_str; i++)
        {
                // Check all the parenthesis characters to check if this is a parenthesis character.
                for(int i_par = 0; i_par < n_pars; i_par++)
                {
                        // Check first for pairs that close a frame.
                        if(right_pars[i_par] == str[i])
                        {
                                // Make sure that this does close the last left parenthesis character.
                                if((l_left_par_chars < 1) ||
                                        (l_left_par_chars >= 1 && left_par_chars[l_left_par_chars-1] != left_pars[i_par]))
                                {
                                        // This character may actually be a left parenthesis, whose left and right characters are the same.
                                        if(left_pars[i_par] == right_pars[i_par])
                                        {
                                                // Interpret this as the left parenthesis characters.
                                                left_par_chars[l_left_par_chars] = str[i];
                                                l_left_par_chars++;
                                        }
                                        else
                                        {
                                                printf("Knotted parenthesis: Right parenthesis character %c is closing left parenthesis character %c.\n", right_pars[i_par], left_par_chars[l_left_par_chars-1]);
                                                return(false);
                                        }
                                }
                                else
                                {
                                        // Note that for the pairs of parenthesis characters whose left and right elements are the same, nesting is not possible.

                                        // Remove this left parenthesis character since its matching right parenthesis is found.
                                        left_par_chars[l_left_par_chars-1] = 0;
                                        l_left_par_chars--;
                                }
                        }
                        else if(left_pars[i_par] == str[i])
                        {
                                // Add this left parenthesis character to the end of list of left parenthesis characters encountered so far.
                                left_par_chars[l_left_par_chars] = str[i];
                                l_left_par_chars++;
                        }
                } // i_par loop.

                // Update par_cnt if this is a parenthesis character, depending on it being left or right parenthesis character.

                // If the last parenthesis character was a left one, the new character cannot be a right one and vice versa. (No knotted pairs.)
        } // i loop

        if(l_left_par_chars == 0)
        {
                return(true);
        }
        else
        {
                return(false);
        }
}

void t_string::clean_string_list(vector<char*>* string_list)
{
	for(int i = 0; i < (int)string_list->size(); i++)
	{
		delete [] string_list->at(i);
	}

	string_list->clear();

	delete string_list;
}

void t_string::clean_tokens(t_string_tokens* tokens)
{
	for(int i = 0; i < (int)tokens->size(); i++)
	{
		delete(tokens->at(i));
	}

	tokens->clear();
	delete(tokens);
}

char* t_string::remove_beginning_spaces_tabs(char* string)
{
	int i_last_space = 0;

	int l_str = t_string::string_length(string);

	for(int i = 0; i < l_str; i++)
	{
		if(string[i] != ' ' &&
			string[i] != '\t')
		{
			i_last_space = i;
			break;
		}
	}

	char* reloc_copy = new char[l_str + 2];

	// Relocate.
	copy(reloc_copy, &string[i_last_space]);
	
	return(reloc_copy);
}

char* t_string::remove_beginning_spaces(char* string)
{
	int i_last_space = 0;

	int l_str = t_string::string_length(string);

	for(int i = 0; i < l_str; i++)
	{
		if(string[i] != ' ')
		{
			i_last_space = i;
			break;
		}
	}

	char* reloc_copy = new char[l_str + 2];

	// Relocate.
	copy(reloc_copy, &string[i_last_space]);
	
	return(reloc_copy);
}

void t_string::remove_beginning_spaces()
{
	int i_last_space = 0;

	for(int i = 0; i < this->length(); i++)
	{
		if(this->x(i) != ' ')
		{
			i_last_space = i;
			break;
		}
	}

	char* reloc_copy = (char*)malloc(sizeof(char*) * (this->length() + 2));

	copy(reloc_copy, this->str());

	// Relocate.
	copy(this->str(), &reloc_copy[i_last_space]);
	free(reloc_copy);
}

t_string_tokens* t_string::tokenize_by_str(char* string, const char* delimiter_string)
{
	t_string_tokens* token_vector = new vector<t_string*>();

	// Current token is nothing at the beginning.
	t_string* current_token = new t_string();

	// Go over the obj_string.
	int str_len = t_string::string_length(string);
	for(int i_str = 0; i_str < str_len; i_str++)
	{
		// This is set to true if i_str'th character is a delimiter.
		bool is_delimiter = false;

		int i_str_search = i_str;

		// Check if current character in the string is a delimiter character.
		for(int i_del = 0; i_del < string_length(delimiter_string); i_del++)
		{
			if(i_str_search == str_len || 
				string[i_str_search] != delimiter_string[i_del])
			{
				is_delimiter = false;
				break;
			}				

			// If the current substring matches the delimiter string, set delimiter to true.
			if(i_del == string_length(delimiter_string) - 1)
			{
				// Push current token to vector (if it is appropriate for pushing.) and start a new token.
				// Push the previous token. Cannot push the next token since there is no way of knowing that
				// it is a valid token beforehand.
				if(current_token != NULL && current_token->length() != 0)
				{
					// Push this token to token vector.
					token_vector->push_back(current_token);

					// Do not free current token now, allocate a new one.
					current_token = new t_string();
				}
				else // The token is not a good token (i.e., an empty token), free its memory and 
				{
					current_token->empty();
				}

				is_delimiter = true;

				// Note that the last char that is at i_str_search is the last char
				// of delimiter string but that will be jumped over by for loop.
				// This is the same in tokenize_by_chars function.
				i_str = i_str_search;
				break;
			}

			// Increment substring counter.
			i_str_search++;

		} // delimiter char loop.

		// If this character is not a delimiter, concatenate the character to the current_token.
		if(!is_delimiter)
		{
			current_token->concat_char(string[i_str]);
		}

	} // string char loop.

	// The last token is added to the tokens vector after whole loop is finished.
	if(current_token != NULL && current_token->length() != 0)
	{
		// Push this token to token vector.
		token_vector->push_back(current_token);
		//printf("pushing %s\n", current_token->obj_string);
	}
	else
	{
		delete(current_token);
	}

	return(token_vector);
}


// Tokenizer: Tokenizes with respect to the characters in the delimiter list, which is a null terminated list.
t_string_tokens* t_string::tokenize_by_str(const char* delimiter_string)
{
	t_string_tokens* token_vector = new vector<t_string*>();

	// Current token is nothing at the beginning.
	t_string* current_token = new t_string();

	// Go over the obj_string.
	int str_len = this->length();
	for(int i_str = 0; i_str < str_len; i_str++)
	{
		//printf("Processing: %c\n", this->x(i_str));
		// This is set to true if i_str'th character is a delimiter.
		bool is_delimiter = false;

		int i_str_search = i_str;
		// Check if current character in the string is a delimiter character.
		for(int i_del = 0; i_del < string_length(delimiter_string); i_del++)
		{
			if(i_str_search == this->length() || this->x(i_str_search) != delimiter_string[i_del])
			{
				is_delimiter = false;
				break;
			}				

			// If the current substring matches the delimiter string, set delimiter to true.
			if(i_del == string_length(delimiter_string) - 1)
			{
				// Push current token to vector (if it is appropriate for pushing.) and start a new token.
				// Push the previous token. Cannot push the next token since there is no way of knowing that
				// it is a valid token beforehand.
				if(current_token != NULL && current_token->length() != 0)
				{
					// Push this token to token vector.
					//printf("pushing %s\n", current_token->obj_string);
					token_vector->push_back(current_token);

					// Do not free current token now, allocate a new one.
					current_token = new t_string();
				}
				else // The token is not a good token (i.e., an empty token), free its memory and 
				{
					current_token->empty();
				}

				is_delimiter = true;

				// Note that the last char that is at i_str_search is the last char
				// of delimiter string but that will be jumped over by for loop.
				// This is the same in tokenize_by_chars function.
				i_str = i_str_search;
				break;
			}

			// Increment substring counter.
			i_str_search++;

		} // delimiter char loop.

		// If this character is not a delimiter, concatenate the character to the current_token.
		if(!is_delimiter)
		{
			current_token->concat_char(this->x(i_str));
		}

	} // string char loop.

	// The last token is added to the tokens vector after whole loop is finished.
	if(current_token != NULL && current_token->length() != 0)
	{
		// Push this token to token vector.
		token_vector->push_back(current_token);
		//printf("pushing %s\n", current_token->obj_string);
	}
	else
	{
		delete(current_token);
	}

	return(token_vector);
}

// Variable argumented sprintf: Print the string of this object using fmt_string as formatter string.
void t_string::sprintf(char* fmt_string, ...)
{
	va_list pars;
	va_start(pars, fmt_string);

	// Initialize string as having nothing.
	this->copy("");

	// Go over the string: Every time a formatting string is found by %.., it is processed.
	for(int i_str = 0; i_str < string_length(fmt_string); i_str++)
	{
		if(fmt_string[i_str] == '%')
		{
			// Check the formatter character.
			i_str++;

			if(fmt_string[i_str] == 'd')
			{
				int _num = va_arg(pars,int);
				this->concat_int(_num);
			}
			else if(fmt_string[i_str] == 'c')
			{
				int _char = va_arg(pars,int);
				this->concat_char(_char);
			}
			else if(fmt_string[i_str] == 's')
			{
				char* _str = (char*)va_arg(pars, void*);
				this->concat_string(_str);
			}
			else if(fmt_string[i_str] == '%')
			{
				this->concat_char(fmt_string[i_str]);
			}

			// The for loop jumps over formatter character.
			//i_str++;
		}
		else // This is a regular character, just copy it.
		{
			this->concat_char(fmt_string[i_str]);
		}
	}

	va_end(pars);
}

void t_string::concat_char(char _char)
{
	//fprintf(stderr, "Concatting %c\n", _char);

	// A memory check.
	int str_len = this->length();
	while(this->string_buffer->l_buffer_mem <= (str_len + 10))
	{		
		char* temp_buf = this->string_buffer->string_buff;
		
		this->string_buffer->l_buffer_mem *= 2;
		this->string_buffer->i_buff = 0;
		this->string_buffer->string_buff = new char[this->string_buffer->l_buffer_mem];		

		// Following uses concat char, again, but at this point, it is for sure that the string has enough memory.
		for(int i = 0; i < str_len; i++)
		{
			this->concat_char(temp_buf[i]);
		}

		delete [] temp_buf;
	}

	int unconcat_length = str_len;
	this->string_buffer->string_buff[unconcat_length] = _char;
	this->string_buffer->i_buff = unconcat_length+1;

	// This is vert important. Finish the string!
	this->string_buffer->string_buff[unconcat_length + 1] = 0;
}

void t_string::concat_string(char* string)
{
	//// Add all the chars in the string.
	//int str_len = string_length(string);
	
	//fprintf(stderr, "Concatting: %s\n", string);
	int i = 0;
	while(string[i] != 0)
	{
		this->concat_char(string[i]);
		i++;
	}
}

void t_string::concat_string(t_string* string)
{
	this->concat_string(string->string_buffer->string_buff);
}

void t_string::concat_int(int i_num)
{
	t_string* num_str = num2str(i_num, 10);
	this->concat_string(num_str);
	delete(num_str);
}

void t_string::concat_float(double f_num)
{
}

// Take the reverse of the string in this object's string buffer.
void t_string::revert()
{
	t_string* temp_buf = new t_string(this->str());
	int str_len = temp_buf->length();
	for(int i = 0; i < str_len; i++)
	{
		this->x(i) = temp_buf->x(str_len - i - 1);
	}

	delete(temp_buf);
}

char* t_string::revert(char* string)
{
	char* reverted_string = t_string::copy_me_str(string);

	//t_string* temp_buf = new t_string(this->str());
	int str_len = t_string::string_length(string);
	for(int i = 0; i < str_len; i++)
	{
		reverted_string[i] = string[str_len - i - 1];
	}

	return(reverted_string);
}

bool t_string::compare(t_string* string)
{
	return(this->compare(string->str()));
}

bool t_string::compare_strings(t_string* str1, t_string* str2)
{
	return(compare_strings(str1->str(), str2->str())); 
}

bool t_string::compare_prefixes(const char* str1, const char* str2)
{
	int i = 0;
	while(str1[i] && str2[i])
	{
		if(str1[i] != str2[i])
		{
			return(false);
		} // comparison.

		i++;
	} // i loop.

	return(true);
}

bool t_string::compare_strings(const char* str1, const char* str2)
{
	if(str1 == NULL || str2 == NULL)
	{
		return(false);
	}

	int str_len1 = string_length(str1);
	int str_len2 = string_length(str2);
	if(str_len1 != str_len2)
	{
		return(false);
	}

	for(int i = 0; i < str_len1; i++)
	{
		if(str1[i] != str2[i])
		{
			return(false);
		}
	}

	return(true);
}

bool t_string::compare_strings_ci(t_string* str1, t_string* str2)
{
	return(compare_strings_ci(str1->str(), str2->str()));
}

bool t_string::compare_substrings_ci(const char* str, const char* sub_str, int& i_match)
{
	// Go over all the substring of str and search for sub_str.
	int i_max = (t_string::string_length(str) - t_string::string_length(sub_str) + 1);

	int l_sub = t_string::string_length(sub_str);
	for(int i = 0; i < i_max; i++)
	{
		char* cur_substr = t_string::substring(str, i, i+l_sub-1);

		if(cur_substr == NULL)
		{
			return(false);
		}

		if(t_string::compare_strings_ci(cur_substr, sub_str))
		{
			i_match = i;
			return(true);
		}

		delete [] cur_substr;

	} // i loop.

	return(false);
}

bool t_string::compare_strings_ci(const char* str1, const char* str2)
{
        int str_len1 = string_length(str1);
        int str_len2 = string_length(str2);

        if(str_len1 != str_len2)
        {
                return(false);
        }

        for(int i = 0; i < str_len1; i++)
        {
                if(toupper(str1[i]) != toupper(str2[i]))
                {
                        return(false);
                }
        }

        return(true);
}

bool t_string::compare_ci(const char* string)
{
        return(compare_strings_ci(this->str(), string));
}

bool t_string::compare_ci(t_string* string)
{
        return(compare_strings_ci(this->str(), string->str()));
}


bool t_string::compare(const char* string)
{
	return(compare_strings(this->str(), string));
}

bool t_string::ends_with(const char* string)
{
	if(this->length() < string_length(string))
	{
		return(false);
	}

	// Go over all the characters of the sub-string and check if they match to the end of the string.
	int i_str = this->length() - 1;
	unsigned int i_sub = string_length(string) - 1;
	while(1)
	{
		if(string[i_sub] != this->x(i_str))
		{
			return(false);
		}

		if(i_sub == 0)
		{
			return(true);
		}
		else
		{
			i_sub--;
		}

		i_str--;		
	} // i loop.

	return(true);
}

bool t_string::starts_with(const char* string)
{
	if(this->length() < string_length(string))
	{
		return(false);
	}

	for(int i_str = 0; i_str < string_length(string); i_str++)
	{
		if(this->x(i_str) != string[i_str])
		{
			return(false);
		}
	}	

	return(true);
}

bool t_string::prefix_matches(const char* str1, const char* str2, bool case_sensitive)
{
	int l_match = 0;
	int l1 = string_length(str1);
	int l2 = string_length(str2);
	if(l1 < l2)
	{
		l_match = l1;
	}
	else
	{
		l_match = l2;
	}

	for(int i_str = 0; i_str < l_match; i_str++)
	{
		if(case_sensitive)
		{
			if(str1[i_str] != str2[i_str])
			{
				return(false);
			}
		}
		else
		{
			if(toupper(str1[i_str]) != toupper(str2[i_str]))
			{
				return(false);
			}
		}
	}	

	return(true);
}

bool t_string::ends_with(const char* full_string, const char* sub_str)
{
	if(string_length(full_string) < string_length(sub_str))
	{
		return(false);
	}

	// Go over all the characters of the sub-string and check if they match to the end of the string.
	int i_str = string_length(full_string) - 1;
	unsigned int i_sub = string_length(sub_str) - 1;
	while(1)
	{
		if(sub_str[i_sub] != full_string[i_str])
		{
			return(false);
		}

		if(i_sub == 0)
		{
			return(true);
		}
		else
		{
			i_sub--;
		}

		i_str--;		
	} // i loop.

	return(true);
}

bool t_string::starts_with(const char* full_string, const char* sub_str)
{
	if(string_length(full_string) < string_length(sub_str))
	{
		return(false);
	}

	for(int i_str = 0; i_str < string_length(sub_str); i_str++)
	{
		if(full_string[i_str] != sub_str[i_str])
		{
			return(false);
		}
	}	

	return(true);
}

bool t_string::starts_with(t_string* string)
{
	return(this->starts_with(string->str()));
}


t_string* t_string::num2str(int num, int base)
{
	t_string* num_str = new t_string();
	int divident = num;
	int residual = divident % base;

	do
	{
		num_str->concat_char((char)(residual + 48)); // Convert to ascii value.		
		divident = divident / base; // Shift right.
		residual = divident % base; // Compute the last digit.

		if(residual > 9)
		{
			printf("The residual greater than 9!\n");
		}
	}
	while(divident != 0);

	// Reverse the num string.
	num_str->revert();

	return(num_str);
}

int t_string::str2num(const char* num_str, int base)
{
    int num = 0;
    int extended_base = 1;
	//for(int i = (int)string_length(num_str)-1; i >= 0; i--)
	unsigned int i = string_length(num_str)-1;
	while(1)
    {
        int cur_digit = (int)num_str[i];
		if(cur_digit >= '0' && cur_digit <= '9')
					num += extended_base * (cur_digit - (int)'0');
		else if(cur_digit >= 'A' && cur_digit <= 'F')
				num += extended_base * (cur_digit - (int)'A' + 10);
		else if(cur_digit >= 'a' && cur_digit <= 'f')
			num += extended_base * (cur_digit - (int)'a' + 10);
		else
		{
			printf("Could not resolve character as number in %s for base %d\n", num_str, base);
			exit(0);
		}

		if(i == 0)
		{
			break;
		}
		else
		{
			i--;
		}

        //printf("%d = %d * (%d - %d)\n", num, extended_base, cur_digit, 48);
        extended_base *= base;
    }

    return(num);
}

void t_string::to_upper()
{
	t_string::to_upper(this->str());
}

void t_string::to_upper(char* string)
{
	//printf("%s->", string);
	//int upper_diff = (int)'A' - (int)'a';

	int str_len = string_length(string); 
	for(int i = 0; i < str_len; i++)
	{
		if(!(string[i] <= 'Z' &&  string[i] >= 'A') && !(string[i] <= 'z' && string[i] >= 'a'))
	 	{
		}
		else if(string[i] <= 'z' && string[i] >= 'a')
		{
			string[i] = string[i] + ((int)'A' - (int)'a');
		}
		else
		{
		}
	}
	//printf("%s\n", string);

}	
int t_string::str2num(t_string* num_str, int base)
{
	return(str2num(num_str->str(), base));
}

