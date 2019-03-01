#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "ansi_string.h"

#ifdef _WIN32
	#include <Windows.h>
#endif

#ifdef unix
	#include <sys/types.h>
	#include <sys/stat.h>
	#include <unistd.h>
#endif

using namespace std;

bool __DUMP_UTIL_MESSAGES__ = false;

char* get_file_extension(char* fp)
{
	t_string_tokens* tokens = t_string::tokenize_by_chars(fp, "/\\.");
	if(tokens->size() > 0)
	{
		char* fn = t_string::copy_me_str(tokens->back()->str());
		t_string::clean_tokens(tokens);
		return(fn);
	}
	else
	{
		return(NULL);
	}
}

bool does_file_exist(char* path)
{
#ifdef unix
	struct stat stat_var;
	int ret = stat(path, &stat_var);

	if(ret == 0)
	{
		return(true);
	}
	else
	{
		return(false);
	}
#endif 

#ifdef _WIN32
	fprintf(stderr, "no stat() implemented in Windows.\n");
	return(false);
#endif
}

void create_directory(char* directory_path)
{
#ifdef _WIN32
	CreateDirectory((LPCWSTR)directory_path, 0);
#endif

#ifdef unix
	mkdir(directory_path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif
}

void get_absolute_fp(char* fp, char* abs_fp)
{
	char cwd[1000];
	get_current_working_directory(cwd);

	if(fp[0] == '/')
	{
		strcpy(abs_fp, fp);
	}
	else
	{
		sprintf(abs_fp, "%s/%s", cwd, fp);
	}
}

void get_current_working_directory(char* buffer)
{
#ifdef _WIN32
	#error "Cannot compile getcwd.";
#endif

#ifdef unix
	char wd_buff[1000];
	getcwd(wd_buff, 1000);
	strcpy(buffer, wd_buff); 
#endif
}

void remove_directory(char* directory_path)
{
}

char* get_file_name(char* fp)
{
	t_string_tokens* tokens = t_string::tokenize_by_chars(fp, "/\\");
	if(tokens->size() > 0)
	{
		char* fn = t_string::copy_me_str(tokens->back()->str());
		t_string::clean_tokens(tokens);
		return(fn);
	}
	else
	{
		return(NULL);
	}
}

char* get_directory_per_fp(char* fp)
{
	t_string_tokens* fp_tokens = t_string::tokenize_by_chars(fp, "\\/");

	if(fp_tokens->size() == 1)
	{
		t_string::clean_tokens(fp_tokens);
		return(NULL);
	}
	else
	{
		char* directory = new char[strlen(fp) + 2];
		directory[0] = 0;

		for(int i = 0; i < (int)fp_tokens->size() - 1; i++)
		{
			strcat(directory, fp_tokens->at(i)->str());
		} // i loop.

		t_string::clean_tokens(fp_tokens);

		return(directory);
	}
}

vector<char*>* load_directory_files(char* root_dir, char* extension)
{
        char ls_cmd[1000];

		if(extension != NULL)
		{
			sprintf(ls_cmd, "ls -d %s/*.%s | xargs -Ifiles basename files > cmd_op.txt", root_dir, extension);
			system(ls_cmd);
		}
		else
		{
			sprintf(ls_cmd, "ls -d %s/* | xargs -Ifiles basename files > cmd_op.txt", root_dir);
			system(ls_cmd);
		}

        vector<char*>* dir_files = new vector<char*>();

        FILE* f_cmd_op = fopen("cmd_op.txt", "r");
        char current_fn[1000];
        while(fscanf(f_cmd_op, "%s", current_fn) == 1)
        {
                char* new_fn = new char[strlen(current_fn) + 3];
                strcpy(new_fn, current_fn);
                dir_files->push_back(new_fn);
        }
        fclose(f_cmd_op);

		// Erase the file.
        sprintf(ls_cmd, "rm -f cmd_op.txt");
        system(ls_cmd);

        return(dir_files);
}

long int get_file_size(char* fp)
{
	FILE* f = NULL;

	// If the file is "stdin", load the standard input.
	if(strcmp(fp, "stdin") == 0)
	{
		fprintf(stderr, "Cannot buffer stdin.\n");
		return(-1);
	}
	else
	{
		// Open the file if it is not the stdin.
		f = fopen(fp, "rb");
		if (f == NULL) 
		{ 
			fprintf(stderr, "Could not open %s\n", fp);
			return(-1);
		} 
	}

	fseek(f, 0, SEEK_END);
	long int file_size = ftell(f);
	fprintf(stderr, "Loading %s of %ld bytes.\n", fp, file_size);

	return(file_size);
}

// Loads a file into memory.
t_file_buffer* load_file(char* fp)
{
	FILE* f = NULL;

	// If the file is "stdin", load the standard input.
	if(strcmp(fp, "stdin") == 0)
	{
		fprintf(stderr, "Cannot buffer stdin.\n");
		return(NULL);
	}
	else
	{
		// Open the file if it is not the stdin.
		f = fopen(fp, "rb");
		if (f == NULL) 
		{ 
			fprintf(stderr, "Could not open %s\n", fp);
			return(NULL);
		} 
	}

	fseek(f, 0, SEEK_END);
	long int file_size = ftell(f);
if(__DUMP_UTIL_MESSAGES__)
{
	fprintf(stderr, "Loading %s of %ld bytes.\n", fp, file_size);
}

	fseek(f, 0, SEEK_SET);
	char* file_buff = new char[file_size + 1];

if(__DUMP_UTIL_MESSAGES__)
{
	fprintf(stderr, "Allocated buffer memory.\n");
}

	if (file_size != (long int)fread(file_buff, sizeof(char), file_size, f)) 
	{ 
		delete [] file_buff;
		return NULL;
	}

if(__DUMP_UTIL_MESSAGES__)
{
	fprintf(stderr, "Read successfully.\n");
}

	fclose(f);
	file_buff[file_size] = 0;

	t_file_buffer* file_buffer = new t_file_buffer();
	file_buffer->file_buffer = file_buff;
	file_buffer->f_ptr = 0;
	file_buffer->l_file = file_size;

	return file_buffer;
}

void unload_file(t_file_buffer* file_buff)
{
	delete [] file_buff->file_buffer;
	delete file_buff;
}

int get_next_string(FILE* file, char* buff, int buff_size)
{
	memset(buff, 0, buff_size);

	int i_s = 0;
	int ret = 0;
	ret = getc(file);
	while(ret != ' ' && ret != EOF && ret != '\n' && ret != '\t')
	{		
		buff[i_s] = ret;
		i_s++;
		ret = getc(file);
	} // file reading loop.

	return(ret);
}

vector<char*>* buffer_file(char* fp)
{
	vector<char*>* file_lines = new vector<char*>();
	//printf("Buffering %s.\n", fp);

	int n_lines = 0;
	FILE* f = NULL;

	if(strcmp(fp, "stdin") == 0)
	{
		f = stdin;	
	}
	else
	{
		f = fopen(fp, "r");
	}

	if(f == NULL)
	{
		//fprintf(stderr, "Could not open %s\n", fp);
		return(NULL);
	}
	while(1)
	{
		char* new_line = getline(f);
		if(new_line == NULL)
		{
			break;
		}
		file_lines->push_back(new_line);

		n_lines++;

if(__DUMP_UTIL_MESSAGES__)
{
		if((n_lines % 10000) == 0)
		{
			printf("Loaded %d lines.               \r", n_lines);
		}
}
	} // file buffering loop.

if(__DUMP_UTIL_MESSAGES__)
{
	fprintf(stderr, "\n");
}

	fclose(f);

	return(file_lines);
}

bool check_file(char* fp)
{
	FILE* f_temp = fopen(fp, "r");
	if(f_temp == NULL)
	{		
		return(false);
	}

	fclose(f_temp);
	return(true);
}

// Check for valid CR LF's depending on OS,
// Should be run for all ASCII input files.
void validate_file(char* fp)
{
#ifdef _WIN32
	char cur_char;
	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		if(cur_char == CR)
		{
			if(fread(&cur_char, 1, 1, f_ip_bin) == 1)
			{
				if(cur_char != LF)
				{
					// Just a warning here.
					printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
					//exit(0);
				}
			}
			else
			{
				// Just a warning here.
				printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
				//exit(0);
			}
		}
		else if(cur_char == LF) // If there is an immediate LF before seeing a CR, this is a linux file.
		{
			// Just a warning here.
			printf("%s is not compatible with dos ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}

	}
	fclose(f_ip_bin);
#endif

#ifdef __unix__
	char cur_char;
	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		// Linux files do not contain CR's.
		// They only contain LF's.
		if(cur_char == CR)
		{
			// Just a warning here.
			printf("%s is not compatible with Linux ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}
	}
	fclose(f_ip_bin);
#endif

#ifdef __APPLE__
	char cur_char;

	// Open file in binary.
	FILE* f_ip_bin = open_f(fp, "rb");
	while(fread(&cur_char, 1, 1, f_ip_bin) == 1)
	{
		// MAC files do not contain LF's.
		if(cur_char == LF)
		{
			// Just a warning here.
			printf("%s is not compatible with Linux ascii files. CR+LF problem at %s(%d).\n", fp, __FILE__, __LINE__);
			//exit(0);
		}
	}
	fclose(f_ip_bin);
#endif

}


FILE* open_f(const char* fp, const char* mode)
{
	if(fp == NULL || mode == NULL)
	{
		printf("Invalid arguments to open_f for %s.\n", fp);
		exit(0);
	}

	FILE* f = fopen(fp, mode);

	if(f == NULL)
	{
		if(mode[0] == 'r')
		{
			printf("Could not open %s for reading.\n", fp);
			exit(0);
		}
		else if(mode[0] == 'w')
		{
			printf("Could not open %s for writing.\n", fp);
			exit(0);
		}
		else
		{
			printf("Could not open %s for requested operation.\n", fp);
			exit(0);
		}
	}

	return(f);
}

int get_n_non_empty_lines(char* fp)
{
	int n_lines = 0;
	FILE* f = open_f(fp, "r");
	while(1)
	{
		char* cur_line = getline(f);
		if(cur_line == NULL)
		{
			break;
		}

		if(!line_empty(cur_line))
		{
			n_lines++;
		}

		delete [] cur_line;
	}
	fclose(f);

	return(n_lines);
}

bool line_empty(char* line)
{
	int l = strlen(line);
	for(int i = 0; i < l; i++)
	{
		if(line[i] != ' ' &&
			line[i] != '\t' &&
			line[i] != '\n')
		{
			return(false);
		}
	}

	return(true);
}

char* getline_per_file_buffer(t_file_buffer* file_buffer)
{
	//t_string* line_str = new t_string();	
	int i_buff = 0;
	int l_buff = 500;
	char* cur_buff = new char[l_buff];
	memset(cur_buff, 0, l_buff);

	char ret = 0;
	while(1)
	{
		//ret = getc(file);
		char cur_char;
		ret = get_next_char_per_file_buffer(file_buffer, cur_char);

		if(ret == false)
		{
			break;
		}
		
#ifdef _WIN32
		if(cur_char == CR)
		{
			// Make sure that the CR is followed by the correct character on the current system.
			ret = get_next_char_per_file_buffer(file_buffer, cur_char);

			if(ret == false)
			{
				fprintf(stderr, "Could not read a character after CR.\n");
				exit(0);
			}
			else if(cur_char != LF)
			{
				fprintf(stderr, "CR is not followed by LF.\n");
				exit(0);
			}

			break;
		}
#elif defined(__unix__)
		// There is nothing to check in UNIX since a carriage return is a new line.
		if(cur_char == CR)
		{
			fprintf(stderr, "Encountered CR character in __unix__.\n");
			exit(0);
		}
		else if(cur_char == LF)
		{
			// Unix has LF as new line. 
			break;
		}
#elif defined(__APPLE__)
		if(cur_char == LF)
		{
			fprintf(stderr, "Encountered CR character in __unix__.\n");
			exit(0);
		}
		else if(cur_char == CR)
		{
			// Unix has LF as new line. 
			break;
		}		
#else
#error "Neither _WIN32 nor __unix__ nor __APPLE__ is not defined.\n";
#endif

		// Update the buffer if it became too long.
		if(i_buff > (l_buff - 5))
		{
			l_buff *= 2;
			char* new_buff = new char[l_buff];
			memset(new_buff, 0, l_buff);
			strcpy(new_buff, cur_buff);
			delete [] cur_buff;
			cur_buff = new_buff;
		}

		cur_buff[i_buff] = cur_char;
		i_buff++;
	} // file reading loop.

	// If the end-of-file is encountered and there was nothing read, return NULL to indicate the there is nothing left to read and EOF is reached.
	if(ret == false && 
		//line_str->length() == 0)
		i_buff == 0)
	{
		//delete(line_str);
		delete [] cur_buff;
		return(NULL);
	}

	return(cur_buff);	
}

bool get_next_char_per_file_buffer(t_file_buffer* file_buffer, char& char_val)
{
	if(file_buffer->f_ptr < file_buffer->l_file)
	{
		char_val = file_buffer->file_buffer[file_buffer->f_ptr];
		file_buffer->f_ptr++;
		return(true);
	}
	else
	{
		char_val = 0;
		return(false);
	}
}

int get_n_lines(FILE* file)
{
	int n_lines = 0;

	char ret = 0;
	int n_chars_in_cur_line = 0;
	while(1)
	{
		ret = getc(file);

		if(ret == EOF)
		{
			break;
		}
		else if(ret == '\n')
		{
			if(n_chars_in_cur_line)
				n_lines++;

			n_chars_in_cur_line = 0;
		}
		else
		{
			n_chars_in_cur_line = 1;
		}
	} // file reading loop.

	return(n_lines);
}

bool get_next_token_per_file_till_newline(FILE* file, char* buffer, int l_buffer, char* delim_char, bool& new_line_check, bool& eof_check)
{
	int i_buff = 0;
	char ret = 0;
	new_line_check = false;
	eof_check = false;
	while(1)
	{
		ret = getc(file);

		// If we found the end-of-file, break, then check whether there is anything in the buffer.
		if(ret == EOF)
		{
			eof_check = true;
			break;
		}

		// If we found a new line, break. The file pointer points to the new line's beginning.
		if(ret == '\n')
		{
			new_line_check = true;
			break;
		}

		// The current character is not a new line, is it the delimiter?
		if(ret == delim_char[0])
		{
			break;
		}
		else if(i_buff < l_buffer)
		{
			// Copy the read character to the buffer, if the buffer did not overflow, otherwise does not copy but keeps reading the file till the next token is found.
			buffer[i_buff] = ret;
			i_buff++;
		}
	} // file reading loop.

	if(i_buff < l_buffer)
	{
		buffer[i_buff] = 0;
	}
	else
	{
		buffer[l_buffer-1] = 0;
	}
	
	if(i_buff == 0)
	{
		return(false);
	}
	else
	{
		return(true);
	}
}

char* getline(FILE* file)
{
	vector<char>* line_vec = new vector<char>();

	char ret = 0;
	while(1)
	{
		ret = getc(file);

		if(ret == EOF)
		{
			break;
		}
		else if(ret == '\n')
		{
			break;
		}

		// Add this character to the vector of characters.
		line_vec->push_back(ret);
	} // file reading loop.

	// If the end-of-file is encountered and there was nothing read, return NULL to indicate the there is nothing left to read and EOF is reached.
	if(ret == EOF && 
		line_vec->size() == 0)
	{
		delete(line_vec);
		return(NULL);
	}

	char* buff = new char[line_vec->size() + 2];
	memset(buff, 0, sizeof(char) * (line_vec->size() + 2));

	// Copy the vector to buffer.
	for(int i = 0; i < (int)line_vec->size(); i++)
	{
		buff[i] = line_vec->at(i);
	} // i loop.

	delete(line_vec);

	return(buff);
}

char* resolve_data_dir()
{
	// try to resolve the DATAPATH_ENV_VAR.
	char* data_dir_from_env = getenv(DATAPATH_ENV_VAR);

	if(data_dir_from_env != NULL)
	{
		char* data_dir = t_string::copy_me_str(data_dir_from_env);
		return(data_dir);
	}
	else
	{
		char* data_dir = t_string::copy_me_str(LOCAL_DATA_PATH);
		return(data_dir);
	}

	printf("Could not resolve thermodynamics data directory.\n");
	exit(0);
}

char* x_fgets(char* buff, int size, FILE* file)
{
	if(fgets(buff, size, file) == NULL)
	{
		return(NULL);
	}

	if(buff[strlen(buff) - 1] == '\n')
	{
		int i_new_line_char = (int)strlen(buff) - 1;
		buff[i_new_line_char] = 0;
	}

	return(buff);
}

