#ifndef _ANSI_STRING_
#define _ANSI_STRING_

//#include <iostream>
#include <vector>
using namespace::std;

// This is the initial buffer size for string with no initializer.
#define IBS (500)

class t_string;

typedef vector<t_string*> t_string_tokens;

/*
The ansi string library:
Portable string representation for extendible, memory managed string with big capabilities,
This is necessary for not having to implement a new version of this code every time a new application is
implemented. tokenizer and sprintf are two vital functions to implement.
*/
struct t_str_node
{
	char* str;	
	int i;
	void* dat;
};

// This is for fast processing of strings. This stores all the information for fast processing of strings. The string buffer should be replaced in t_string with this struct.
struct t_string_buffer
{
	char* string_buff;
	int l_buffer_mem;
	int i_buff; // The end position of the string in the buffer. This enables returning length of the string in O(1) time.
};

class t_string
{
public:
	t_string(char* string);
	t_string(t_string* string);
	t_string();
	~t_string();

	static void set_byte_buffer(void* buffer, long int n_chars_2_set, char val_2_set);

	static vector<int>* get_string_sorting_idx(vector<char*>* string_list);
	static bool sort_str_nodes(t_str_node* node1, t_str_node* node2);

	static char* copy_me_str(const char* str);

	static int fast_count_non_empty_tokens(char* str, char* delims);

	static bool sort_strings(char* str1, char* str2);
	static bool sort_strings_per_prefix(char* str1, char* str2);
	static int fast_search_string(char* query, vector<char*>* string_list, int i, int j);
	static int fast_search_string_per_prefix(char* query, vector<char*>* string_list, int i, int j);

	//char* obj_string; // null terminated string that contains the string information.
	//unsigned int obj_str_mem_size; // memory size of the string, not the actual size. This number is greater than or equal to length of the string.

	// This is the heart of all processing.
	t_string_buffer* string_buffer;

	static void clean_string_list(vector<char*>* string_list);

	static char* char_2_str(const char cur_char);

	// Tokenize the string of this string object and return a
	// vector of string objects.
	t_string_tokens* tokenize_by_chars(const char* delimiter_list);
	t_string_tokens* tokenize_by_str(const char* delimiter_string);
	static t_string_tokens* tokenize_by_chars(char* string, const char* delimiter_list, bool return_empty_tokens);
	static t_string_tokens* tokenize_by_chars(char* string, const char* delimiter_list);
	static t_string_tokens* tokenize_by_str(char* string, const char* delimiter_string);
	static vector<char*>* copy_tokens_2_strs(t_string_tokens* toks);
	char* substring(int i, int j);

	void copy(const char* string);
	void copy(t_string* string);

	static vector<t_string*>* get_first_n_tokens(char* string, int n_tokens, const char* delimiter_list, int& i_next_char);
	static bool get_next_token(char* string, char* buffer, int l_buffer, const char* delimiter_list, int& i_cur_char);

	static vector<vector<char*>*>* split_string_list(vector<char*>* string_list, int n_batches);
	static vector<vector<void*>*>* split_object_list(vector<void*>* obj_list, int n_batches);
	static vector<vector<int>*>* split_int_list(vector<int>* int_list, int n_batches);
	static vector<void*>* interleave_obj_list(vector<vector<void*>*>* obj_lists_list);
	static vector<void*>* pool_obj_list(vector<vector<void*>*>* obj_lists_list);

	static vector<char*>* get_unique_entries(vector<char*>* str_list);

	// Static string library functions.
	static char* substring(const char* str, int i, int j);
	static void copy(char* dest_string, const char* src_string);
	static int string_length(const char* string);
	static int string_length(t_string* string);
	static t_string* num2str(int num, int base);
	static int str2num(const char* num_str, int base);
	static int str2num(t_string* num_str, int base);
	static bool compare_strings(t_string* str1, t_string* str2);
	static bool compare_prefixes(const char* str1, const char* str2);
    static bool compare_strings(const char* str1, const char* str2);
    static bool compare_strings_ci(t_string* str1, t_string* str2);
    static bool compare_strings_ci(const char* str1, const char* str2);
	static bool compare_substrings_ci(const char* str, const char* sub_str, int& i_match); // Check if the string contains sub_str as a substring in it, without regard to the case of letters.
	static void to_upper(char* string);
	static void clean_tokens(t_string_tokens* tokens);

	static int get_matching_prefix_length(char* str1, char* str2);

	static bool compare_strings_per_total_prefix(char* str1, char* str2);

	static int get_add_i_str(vector<char*>* strs, char* str);
	static int get_i_str(vector<char*>* strs, char* str);
	static int get_i_str_ci(vector<char*>* strs, char* str);
	static int get_i_str(vector<t_string*>* strs, char* str);

	static bool starts_with(const char* full_string, const char* sub_str);
	static bool ends_with(const char* full_string, const char* sub_str);
	static bool prefix_matches(const char* str1, const char* str2, bool case_sensitive);

	static char* get_alphanumeric(char* string);

	static bool is_balanced(char* str, char* left_pars, char* right_pars);

	static void replace_avoid_list(char* str, char start_char, char end_char, char char_to_replace);
	static void replace_avoid_list(char* str, char* avoided_char_list, char char_to_replace);

	static bool is_number(char* string);

	static bool is_empty(char* string);

	static char* get_longest_matching_substring(char* string1, char* string2);

	// Parse the consercutive numbers in the string and return them in a vector.
	static vector<int>* get_integers_in_string(char* string);
	vector<int>* get_integers_in_string();

	void remove_beginning_spaces();
	static char* remove_beginning_spaces(char* string);
	static char* remove_beginning_spaces_tabs(char* string);

	bool is_balanced(char* left_pars, char* right_pars);

	// Concatenation functions. These are the basis of other functions because
	// these do the memory scaling. No other functions should need the memory scaling for the string.
	void concat_char(char _char);
	void concat_string(char* string);
	void concat_string(t_string* string);
	void concat_int(int i_num);
	void concat_float(double f_num);

	// sprintf declaration.
	void sprintf(char* fmt_string, ...);

	bool compare(const char* string);
	bool compare(t_string* string);
	bool compare_ci(t_string* string);
	bool compare_ci(const char* string);
	bool starts_with(const char* string);
	bool starts_with(t_string* string);
	bool ends_with(const char* string);

	char& x(int i);
	char* str();
	int length();
	void empty();
	void revert();	
	static char* revert(char* string);
	void to_upper();
};

#endif // _ANSI_STRING_


