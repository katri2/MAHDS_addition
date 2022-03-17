#include "config.h"


void delete_space_symbols(char* s) {
    int end = 0;
    for (int i = 0; i < strlen(s); i++) {
        if (s[i] != ' ' && s[i] != '\n' && s[i] != '\r' && s[i] != '\t') {
            if (i != end) {
                s[end] = s[i];
            }
            end++;
        }
    }
    s[end] = '\0';
}

void replace_sym (char* s, char a, char b) {
	int i = 0;
	while (s[i]!='\0') {
		if (s[i] == a) {
			s[i] = b;
		}
		++i;
	}
}

void load_config(
	char*           system_settings_path,
	char*           alignment_settings_path,
	size_t          symbol_types,
	size_t          spectre_middle,
	OUT mahds_config* config
) { //alrady got seq
	FILE* config_file = NULL;
	char str[SYMBOLS_IN_CFG_STR], k[SYMBOLS_IN_CFG_STR], v[SYMBOLS_IN_CFG_STR];
	unsigned short got_params = 0;

	if ((config_file = fopen(system_settings_path, "r")) == NULL){
		printf("Can not open system config file\n");
		return;
	}

	while (fgets (str, SYMBOLS_IN_CFG_STR, config_file) != NULL) {
		delete_space_symbols(str);
		replace_sym(str, ':', ' ');
		int got_fields_n = sscanf(str,"%s%s",k,v);
		if (got_fields_n!=2) {
			continue;
		}

		if (strcmp(k,"threads_per_node") == 0) {
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->threads_num = THREADS_NUM;
			} else {
				got_fields_n = sscanf(v,"%ld",&(config->threads_num));
				if (! got_fields_n) {
					printf("Error while parsing 'threads_per_node'\n");
					fclose(config_file);
					return;
				}
			}
		} else if (strcmp(k,"matrix_path") == 0) {
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->matrix_path = malloc((strlen(DEFAULT_MATR_PATH) + 1) * sizeof(char));
				strcpy(config->matrix_path, DEFAULT_MATR_PATH);
			} else {
				config->matrix_path = malloc((strlen(v) + 1) * sizeof(char));
				strcpy(config->matrix_path, v);
			}
		}
	}
	fclose(config_file);
	if (got_params != 2) {
		printf("Error while parsing system configuration: %u fields got\n", got_params);
		return;
	}

	////////////////////////////////////////////////////////////////////////////////////
	got_params = 0;
	char symbol_kinds_got = 0;
	if ((config_file = fopen(alignment_settings_path, "r")) == NULL){
		printf("Can not open alignment config file\n");
		return;
	}

	while (fgets (str, SYMBOLS_IN_CFG_STR, config_file) != NULL) {
		delete_space_symbols(str);
		replace_sym(str, ':', ' ');
		int got_fields_n = sscanf(str,"%s%s",k,v);
		if (got_fields_n!=2) {
			continue;
		}

		if (strcmp(k,"sym_kinds_num") == 0) {   //1
			symbol_kinds_got = 1;
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->size_of_symbols_set = PROT_SYMBOLS;
			} else {
				got_fields_n = sscanf(v,"%ld",&(config->size_of_symbols_set));
				if (! got_fields_n) {
					printf("Error while parsing 'sym_kinds_num'\n");
					fclose(config_file);
					return;
				}
			}
		} else if (strcmp(k,"d_price") == 0) {   //2
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->d_price = D_PRICE;
			} else {
				got_fields_n = sscanf(v,"%lf",&(config->d_price));
				if (! got_fields_n) {
					printf("Error while parsing 'd_price'\n");
					fclose(config_file);
					return;
				}
			}
		} else if (strcmp(k,"e_price") == 0) {   //3
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->e_price = E_PRICE;
			} else {
				got_fields_n = sscanf(v,"%lf",&(config->e_price));
				if (! got_fields_n) {
					printf("Error while parsing 'e_price'\n");
					fclose(config_file);
					return;
				}
			}
		} else if (strcmp(k,"kd") == 0) {   //4
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->K_d = K_D_CONST;
			} else {
				got_fields_n = sscanf(v,"%lf",&(config->K_d));
				if (! got_fields_n) {
					printf("Error while parsing 'kd'\n");
					fclose(config_file);
					return;
				}
			}
		} else if (strcmp(k,"r_mult") == 0) {   //5
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->R_sqr_mult_const = R_MULT_CONST;
			} else {
				got_fields_n = sscanf(v,"%lf",&(config->R_sqr_mult_const));
				if (! got_fields_n) {
					printf("Error while parsing 'r_mult'\n");
					fclose(config_file);
					return;
				}
			}
		} else if (strcmp(k,"borders_width") == 0) {   //6
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->alignment_matrix_borders = BORDERS;
			} else {
				got_fields_n = sscanf(v,"%ld",&(config->alignment_matrix_borders));
				if (! got_fields_n) {
					printf("Error while parsing 'borders_width'\n");
					fclose(config_file);
					return;
				}
			}
		} else if (strcmp(k,"checking_matr_set") == 0) {   //7
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->size_of_checking_matrix_set = CHECKING_SET;
			} else {
				got_fields_n = sscanf(v,"%ld",&(config->size_of_checking_matrix_set));
				if (! got_fields_n) {
					printf("Error while parsing 'checking_matr_set'\n");
					fclose(config_file);
					return;
				}
			}
		} else if (strcmp(k,"spectre_scatter") == 0) {   //8
			++got_params;
			size_t delta;
			if (strcmp(v,"default") == 0) {
				delta = PERIODS_SCATTER;
			} else {
				got_fields_n = sscanf(v,"%ld",&delta);
				if (! got_fields_n) {
					printf("Error while parsing 'spectre_scatter'\n");
					fclose(config_file);
					return;
				}
			}

			if ((long long)spectre_middle-(long long)delta >= 2) {
				config->spectre_min = spectre_middle-delta;
			} else {
				config->spectre_min = 2;
			}

			if (config->spectre_min <= spectre_middle+delta) {
				config->spectre_max = spectre_middle+delta;
			} else {
				config->spectre_max = config->spectre_min;
			}

		} else if (strcmp(k,"alignment_type") == 0) {   //9
			++got_params;
			if (strcmp(v,"default") == 0) {
				config->alignment_type = MAHDS_ALIGNMENT_TYPE_GLOBAL;
			} else if (strcmp(v,"global") == 0) { 
				config->alignment_type = MAHDS_ALIGNMENT_TYPE_GLOBAL;
			} else if (strcmp(v,"local") == 0) { 
				config->alignment_type = MAHDS_ALIGNMENT_TYPE_LOCAL;
			} else if (strcmp(v,"mixed") == 0) { 
				config->alignment_type = MAHDS_ALIGNMENT_TYPE_MIXED;
			} else {
				printf("Error while parsing 'alignment_type'\n");
				fclose(config_file);
				return;
			}
		}
	}
	fclose(config_file);
	if (! (got_params == 9 || (got_params == 8 && !symbol_kinds_got)) ) {
		printf("Error while parsing alignment configuration: %u fields got\n", got_params);
		return;
	}
	if (!symbol_kinds_got) {
		config->size_of_symbols_set = symbol_types;
	}
}
