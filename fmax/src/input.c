#include "input.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
sequences_info* fill_artificial_seq_info (kseq_t* original_set, char with_gaps){ // with_gaps:bool
    sequences_info* res = (sequences_info*)malloc(sizeof(sequences_info));
    res->seq_set = (mahds_seq_set*)malloc(sizeof(mahds_seq_set));
    res->seq_set->seq_number = 0;
    res->seq_set->sequences = NULL;
    res->mapping_info = malloc(sizeof(mapping_info));
    res->mapping_info->original_headers = NULL;
    res->mapping_info->original_set = NULL;
    res->mapping_info->output_gap_symbol = '-';
    char* used_orig_symbols = calloc(sizeof(char), UCHAR_MAX+1);
    mahds_symbol* used_mapped_symbols = calloc(sizeof(mahds_symbol), SYMBOL_MAX+1);

    int l;
    char warning = 0;
    while ((l = kseq_read(original_set)) >= 0) {
        ++res->seq_set->seq_number;
        res->mapping_info->original_headers = (char**)realloc((void*)res->mapping_info->original_headers, res->seq_set->seq_number*sizeof(char*));
        char buf[2048];
        sprintf(buf, ">%s", original_set->name.s);
        if (original_set->comment.l) {
            char bu[2048];
            sprintf(bu, " %s", original_set->comment.s);
            strcat(buf, bu);
        }
        res->mapping_info->original_headers[res->seq_set->seq_number-1] = (char*)malloc(strlen(buf) * sizeof(char));
        strcpy(res->mapping_info->original_headers[res->seq_set->seq_number-1], buf);

        res->seq_set->sequences = (mahds_symbolic_sequence**)realloc((void*)res->seq_set->sequences, res->seq_set->seq_number*sizeof(mahds_symbolic_sequence*));
        res->mapping_info->original_set = (mahds_symbolic_sequence**)realloc((void*)res->mapping_info->original_set, res->seq_set->seq_number*sizeof(mahds_symbolic_sequence*));
        res->seq_set->sequences[res->seq_set->seq_number-1] = (mahds_symbolic_sequence*)malloc(sizeof(mahds_symbolic_sequence));
        res->mapping_info->original_set[res->seq_set->seq_number-1] = (mahds_symbolic_sequence*)malloc(sizeof(mahds_symbolic_sequence));
        res->seq_set->sequences[res->seq_set->seq_number-1]->length = original_set->seq.l;
        res->mapping_info->original_set[res->seq_set->seq_number-1]->length = original_set->seq.l;
        res->seq_set->sequences[res->seq_set->seq_number-1]->symbols = malloc(sizeof(mahds_symbol) * original_set->seq.l);
        res->mapping_info->original_set[res->seq_set->seq_number-1]->symbols = malloc(sizeof(mahds_symbol) * original_set->seq.l);
        for (size_t i = 0; i < original_set->seq.l; ++i) {
            char sym = original_set->seq.s[i];
            res->mapping_info->original_set[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)sym;
            if (sym >= 'A' && sym <= 'Z') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)(sym - 'A');
            } else if (sym >= 'a' && sym <= 'z') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)(sym - 'a');
            // } else if (sym >= '0' && sym <= '9') {
            //     res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)(sym - '0');
            } else if (sym == ' ' || sym == '\t') {
                continue;
            } else if (sym == '-') {
                if (!with_gaps) {
                    warning = 1;
                }
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)SYMBOL_MAX;
            } else if (sym >= '0' && sym <= '9' && with_gaps) {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)sym;
            } else {
                printf("In fill_artificial_seq_info: unexpected symbo! (%c)\n", sym);
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)SYMBOL_MAX;
            }
            used_orig_symbols[sym + abs(CHAR_MIN)] = 1;
            used_mapped_symbols[res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i]] = 1;
        }
	}
    
    ///////
    if (warning) {
        printf("Warning: symbol '-' is present in sequences!\n");
    }
    ///////
    size_t counter = 0;
    for (size_t i = 0; i < (UCHAR_MAX+1); ++i) {
        if (used_orig_symbols[i] == 1) {
            ++counter;
        }
    }
    res->mapping_info->original_count_symbols = counter;
    counter = 0;
    for (size_t i = 0; i < (SYMBOL_MAX+1); ++i) {
        if (used_mapped_symbols[i] == 1) {
            ++counter;
        }
    }
    res->mapping_info->mapped_count_symbols = counter;
    ///////
    if (used_mapped_symbols[SYMBOL_MAX] == 1) {
        --counter;
    }
    for (size_t i = 0; i < res->seq_set->seq_number; ++i){
        res->seq_set->sequences[i]->count_symbols = counter;
    }
    for (size_t i = 0; i < res->seq_set->seq_number; ++i){
        res->mapping_info->original_set[i]->count_symbols = res->mapping_info->original_count_symbols;
    }

    free(used_orig_symbols);
    free(used_mapped_symbols);
    return res;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
sequences_info* fill_DNA_seq_info (kseq_t* original_set, char with_gaps) { // with_gaps:bool
    printf("Function fill_DNA_seq_info is not implemented\n");
    return NULL;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
sequences_info* fill_protein_seq_info (kseq_t* original_set, char with_gaps) { // with_gaps:bool
    sequences_info* res = (sequences_info*)malloc(sizeof(sequences_info));
    res->seq_set = (mahds_seq_set*)malloc(sizeof(mahds_seq_set));
    res->seq_set->seq_number = 0;
    res->seq_set->sequences = NULL;
    res->mapping_info = malloc(sizeof(mapping_info));
    res->mapping_info->original_headers = NULL;
    res->mapping_info->original_set = NULL;
    res->mapping_info->output_gap_symbol = '-';

    char* used_orig_symbols = calloc(sizeof(char), UCHAR_MAX+1);
    mahds_symbol* used_mapped_symbols = calloc(sizeof(mahds_symbol), SYMBOL_MAX+1);

    int l;
    srand(time(NULL));
    while ((l = kseq_read(original_set)) >= 0) {
        ++res->seq_set->seq_number;

        res->mapping_info->original_headers = (char**)realloc((void*)res->mapping_info->original_headers, res->seq_set->seq_number*sizeof(char*));
        char buf[2048];
        sprintf(buf, ">%s", original_set->name.s);
        if (original_set->comment.l) {
            char bu[2048];
            sprintf(bu, " %s", original_set->comment.s);
            strcat(buf, bu);
        }
        res->mapping_info->original_headers[res->seq_set->seq_number-1] = (char*)malloc(strlen(buf) * sizeof(char));
        strcpy(res->mapping_info->original_headers[res->seq_set->seq_number-1], buf);
        //printf("stored: %s\n", res->mapping_info->original_headers[res->seq_set->seq_number-1]);
        //printf("buffer: %s\n____\n", buf);
        res->seq_set->sequences = (mahds_symbolic_sequence**)realloc((void*)res->seq_set->sequences, res->seq_set->seq_number*sizeof(mahds_symbolic_sequence*));
        res->mapping_info->original_set = (mahds_symbolic_sequence**)realloc((void*)res->mapping_info->original_set, res->seq_set->seq_number*sizeof(mahds_symbolic_sequence*));
        res->seq_set->sequences[res->seq_set->seq_number-1] = (mahds_symbolic_sequence*)malloc(sizeof(mahds_symbolic_sequence));
        res->mapping_info->original_set[res->seq_set->seq_number-1] = (mahds_symbolic_sequence*)malloc(sizeof(mahds_symbolic_sequence));
        res->seq_set->sequences[res->seq_set->seq_number-1]->length = original_set->seq.l;
        res->mapping_info->original_set[res->seq_set->seq_number-1]->length = original_set->seq.l;
        res->seq_set->sequences[res->seq_set->seq_number-1]->symbols = malloc(sizeof(mahds_symbol) * original_set->seq.l);
        res->mapping_info->original_set[res->seq_set->seq_number-1]->symbols = malloc(sizeof(mahds_symbol) * original_set->seq.l);
        for (size_t i = 0; i < original_set->seq.l; ++i) {
            char sym = original_set->seq.s[i];
            res->mapping_info->original_set[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)sym;
            if (sym == 'A' || sym == 'a') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)0;
            } else if (sym == 'B' || sym == 'b') {
                if (rand() <= RAND_MAX/2){
                    res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)2;
                } else {
                    res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)11;
                }
            } else if (sym == 'C' || sym == 'c') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)1;
            } else if (sym == 'D' || sym == 'd') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)2;
            } else if (sym == 'E' || sym == 'e') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)3;
            } else if (sym == 'F' || sym == 'f') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)4;
            } else if (sym == 'G' || sym == 'g') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)5;
            } else if (sym == 'H' || sym == 'h') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)6;
            } else if (sym == 'I' || sym == 'i') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)7;
            } else if (sym == 'J' || sym == 'j') {
                if (rand() <= RAND_MAX/2){
                    res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)9;
                } else {
                    res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)7;
                }
            } else if (sym == 'K' || sym == 'k') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)8;
            } else if (sym == 'L' || sym == 'l') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)9;
            } else if (sym == 'M' || sym == 'm') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)10;
            } else if (sym == 'N' || sym == 'n') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)11;
            } else if (sym == 'O' || sym == 'o') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)(rand()%20); // (21th) in fact
            } else if (sym == 'P' || sym == 'p') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)12;
            } else if (sym == 'Q' || sym == 'q') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)13;
            } else if (sym == 'R' || sym == 'r') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)14;
            } else if (sym == 'S' || sym == 's') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)15;
            } else if (sym == 'T' || sym == 't') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)16;
            } else if (sym == 'U' || sym == 'u') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)(rand()%20); // (22th) in fact
            } else if (sym == 'V' || sym == 'v') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)17;
            } else if (sym == 'W' || sym == 'w') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)18;
            } else if (sym == 'Y' || sym == 'y') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)19;
            } else if (sym == 'Z' || sym == 'z') {
                if (rand() <= RAND_MAX/2){
                    res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)3;
                } else {
                    res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)13;
                }
            } else if (sym == 'X' || sym == 'x') {
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)(rand()%20);  // any
            } else if (sym == '*') {
                printf("In fill_protein_seq_info: symbol '*'\n");
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)(rand()%20);
            } else if (sym == '-') {
                if (!with_gaps) {
                    printf("In fill_protein_seq_info: symbol '-'\n");
                }
                res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)SYMBOL_MAX;
            } else {
                if (!with_gaps || sym < '0' || sym > '9'){
                    printf("In fill_protein_seq_info: unexpected symbol: %c\n", sym);
                    res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)SYMBOL_MAX;
                } else {
                    res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i] = (mahds_symbol)sym;
                }
            }
            used_orig_symbols[sym + abs(CHAR_MIN)] = 1;
            used_mapped_symbols[res->seq_set->sequences[res->seq_set->seq_number-1]->symbols[i]] = 1;
        }
	}
    ///////
    size_t counter = 0;
    for (size_t i = 0; i < (UCHAR_MAX+1); ++i) {
        if (used_orig_symbols[i] == 1) {
            ++counter;
        }
    }
    res->mapping_info->original_count_symbols = counter;
    counter = 0;
    for (size_t i = 0; i < (SYMBOL_MAX+1); ++i) {
        if (used_mapped_symbols[i] == 1) {
            ++counter;
        }
    }
    res->mapping_info->mapped_count_symbols = counter;
    ///////

    for (size_t i = 0; i < res->seq_set->seq_number; ++i){
        res->seq_set->sequences[i]->count_symbols = PROT_SYM_COUNT;
    }
    for (size_t i = 0; i < res->seq_set->seq_number; ++i){
        res->mapping_info->original_set[i]->count_symbols = res->mapping_info->original_count_symbols;
    }

    free(used_orig_symbols);
    free(used_mapped_symbols);
    return res;
};


sequences_info* parse_seq_info(char* file_path, input_mapping in_mapping, char with_gaps){ // with_gaps:bool
	gzFile fp = gzopen(file_path, "r");
	kseq_t *seq = kseq_init(fp);
	
    sequences_info* res = NULL;
    if (in_mapping == MAPPING_ARTIFICIAL) {
        res = fill_artificial_seq_info(seq, with_gaps);
    } else if (in_mapping == MAPPING_DNA) {
        res = fill_DNA_seq_info(seq, with_gaps);
    } else if (in_mapping == MAPPING_PROTEIN) {
        res = fill_protein_seq_info(seq, with_gaps);
    } else {
        printf("In parse_seq_info: unexpected input mapping type\n");
    }

	kseq_destroy(seq);
	gzclose(fp);
	return res;
}


void get_data_from_file(
    input_mapping         mapping_type,
    char*                 input_file_path,
    char*                 sys_cfg_path,
    char*                 align_cfg_path,
    OUT mahds_config**      config,
    OUT sequences_info**  seq_inf
){
    (*seq_inf) = parse_seq_info(input_file_path, mapping_type, (char)0);
    double mean_length = calc_mean_length((*seq_inf)->seq_set->seq_number, (*seq_inf)->seq_set->sequences);
    load_config(
        sys_cfg_path,
	    align_cfg_path,
	    (*seq_inf)->seq_set->sequences[0]->count_symbols,
	    mean_length,
	    (*config)
    );
    for (size_t i = 0; i < (*seq_inf)->seq_set->seq_number; ++i) {
        (*seq_inf)->seq_set->sequences[i]->count_symbols = (*config)->size_of_symbols_set;
    }
    if ((*config)->alignment_matrix_borders < mean_length*0.1) {
        printf("Warning: too narrow borders (less then 0.1 * mean length of sequence)\n");
    } else if ((*config)->alignment_matrix_borders > mean_length * 0.5) {
        printf("Note: Wide borders (more then 0.5 * mean length of sequence)\n");
    }
}
