#include "fmax.h"


double get_fmax (char* file_name, input_mapping in_mapping, char* sys_cfg_path, char* align_cfg_path) {
    sequences_info* seq_inf = parse_seq_info(file_name, in_mapping, (char)1); // mapping info empty here
    mahds_seq_set* seq_set = (mahds_seq_set*)malloc(sizeof(mahds_seq_set));
    seq_set->sequences = NULL;
    seq_set->seq_number = 0;
    for (size_t i = 0; i < seq_inf->seq_set->seq_number; ++i) {
        if (!strstr(seq_inf->mapping_info->original_headers[i], "periodicity")) {
            ++seq_set->seq_number;
            seq_set->sequences = (mahds_symbolic_sequence**)realloc((void*)seq_set->sequences, seq_set->seq_number*sizeof(mahds_symbolic_sequence*));
            seq_set->sequences[seq_set->seq_number-1] = copy_mahds_seq(seq_inf->seq_set->sequences[i]);
        }
    }

    ///////////////////////////////////// MAKING ARTIFICIAL ALIGNMENT + config /////////////////////////////////////
    mahds_al_set* S1_set = malloc(sizeof(mahds_al_set));
    S1_set->seq_number = seq_set->seq_number;
    S1_set->sequences = (mahds_alignment_sequence**)malloc(S1_set->seq_number * sizeof(mahds_alignment_sequence*));
    for (size_t i = 0; i < S1_set->seq_number; ++i) {
        S1_set->sequences[i] = (mahds_alignment_sequence*)malloc(sizeof(mahds_alignment_sequence));
        S1_set->sequences[i]->sequence.length = 0;
        S1_set->sequences[i]->sequence.symbols = NULL;
    }
    
    mahds_al_set* S_set = malloc(sizeof(mahds_al_set));
    S_set->seq_number = seq_set->seq_number;
    S_set->sequences = (mahds_alignment_sequence**)malloc(S_set->seq_number * sizeof(mahds_alignment_sequence*));
    for (size_t i = 0; i < S_set->seq_number; ++i) {
        S_set->sequences[i] = malloc(sizeof(mahds_alignment_sequence));
        S_set->sequences[i]->sequence.length = seq_set->sequences[i]->length;
        S_set->sequences[i]->sequence.count_symbols = seq_set->sequences[i]->count_symbols;
        S_set->sequences[i]->sequence.symbols = malloc(S_set->sequences[i]->sequence.length * sizeof(mahds_symbol));
        for (size_t j = 0; j < S_set->sequences[i]->sequence.length; ++j) {
            S_set->sequences[i]->sequence.symbols[j] = seq_set->sequences[i]->symbols[j];
        }
    }
    
    size_t middle = ceil((double)seq_set->seq_number/2.);
    size_t counter;
    mahds_symbol period_sym = 0;
    for (size_t i = 0; i < S_set->sequences[0]->sequence.length; ++i) {
        counter = 0;
        for (size_t s = 0; s < S_set->seq_number; ++s) {
            if (S_set->sequences[s]->sequence.symbols[i] >= S_set->sequences[s]->sequence.count_symbols) { // gaps in column
                ++counter;
            }
        }
        if (counter <= middle) {
            for (size_t s = 0; s < S1_set->seq_number; ++s) {
                ins_symbol(&(S1_set->sequences[s]->sequence), i, period_sym);
            }
            ++period_sym;
        } else {
            for (size_t s = 0; s < S1_set->seq_number; ++s) {
                ins_symbol(&(S1_set->sequences[s]->sequence), i, SYMBOL_MAX);
            }
        }
    }
    
    //////CHECK//////
    if (!period_sym) {
        printf("In get_fmax: period_sym == 0\n");
        exit(EXIT_FAILURE);
    }
    /////////////////

    ////////////////
    mahds_config config;
    load_config(
        sys_cfg_path,
	    align_cfg_path,
	    seq_set->sequences[0]->count_symbols,
	    (size_t)period_sym,
	    &config
    );
    
    for (size_t i = 0; i < seq_set->seq_number; ++i) {
        S_set->sequences[i]->sequence.count_symbols = config.size_of_symbols_set;
        seq_set->sequences[i]->count_symbols = config.size_of_symbols_set;
        S1_set->sequences[i]->sequence.count_symbols = (size_t)period_sym;
    }
    ////////////////

    /////////////// TRANSFORMING seq_set TO SET OF RAW (without gaps) SEQUENCES ////////////////
    for (size_t i = 0; i < seq_set->seq_number; ++i) {
        for (long j = seq_set->sequences[i]->length - 1; j >= 0; --j) {
            if (seq_set->sequences[i]->symbols[j] >= seq_set->sequences[i]->count_symbols) {
                del_symbol(seq_set->sequences[i], j);
            }
        }
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////// CORRECTING ARTIFICIAL ALIGNMENT /////////////////////////////////////
    size_t glob_period_begin = 0;
    size_t glob_period_end = S1_set->sequences[0]->sequence.length; // inx < end
    while (S1_set->sequences[0]->sequence.symbols[glob_period_begin] >= period_sym) {
        ++glob_period_begin;
    }
    while (S1_set->sequences[0]->sequence.symbols[glob_period_end-1] >= period_sym) {
        --glob_period_end;
    }
    
    for (size_t i = 0; i < S_set->seq_number; ++i) {
        if (config.alignment_type != MAHDS_ALIGNMENT_TYPE_GLOBAL) {
            size_t S_real_begin = 0;
            size_t S_real_end = S_set->sequences[i]->sequence.length; // inx < end
            while (S_set->sequences[i]->sequence.symbols[S_real_begin] >= S_set->sequences[i]->sequence.count_symbols) {
                ++S_real_begin;
            }
            while (S_set->sequences[i]->sequence.symbols[S_real_end-1] >= S_set->sequences[i]->sequence.count_symbols){
                --S_real_end;
            }

            if (S_real_begin < glob_period_begin) {
                counter = 0;
                for (size_t k = S_real_begin; k < glob_period_begin; ++k) {
                    if (S_set->sequences[i]->sequence.symbols[k] < S_set->sequences[i]->sequence.count_symbols) {
                        ++counter;
                    }
                }
                S_set->sequences[i]->orig_begin = counter;
            } else {
                counter = 0;
                for (size_t k = glob_period_begin; k < S_real_begin; ++k) {
                    if (S1_set->sequences[i]->sequence.symbols[k] < S1_set->sequences[i]->sequence.count_symbols) {
                        ++counter;
                    }
                }
                S1_set->sequences[i]->orig_begin = counter;
            }

            if (S_real_end < glob_period_end) {
                counter = 0;
                for (size_t k = S_real_end; k < glob_period_end; ++k) {
                    if (S1_set->sequences[i]->sequence.symbols[k] < S1_set->sequences[i]->sequence.count_symbols) {
                        ++counter;
                    }
                }
                S1_set->sequences[i]->orig_end = ((size_t)period_sym - counter);
            } else {
                counter = 0;
                for (size_t k = glob_period_end; k < S_real_end; ++k) {
                    if (S_set->sequences[i]->sequence.symbols[k] < S_set->sequences[i]->sequence.count_symbols) {
                        ++counter;
                    }
                }
                S_set->sequences[i]->orig_end = (seq_set->sequences[i]->length - counter);
            }
            
            // deleting segments
            if (S_real_end < glob_period_end) {
                del_symbols((&S_set->sequences[i]->sequence), S_real_end, S_set->sequences[i]->sequence.length);
                del_symbols((&S1_set->sequences[i]->sequence), S_real_end, S1_set->sequences[i]->sequence.length);
            } else {
                del_symbols((&S_set->sequences[i]->sequence), glob_period_end, S_set->sequences[i]->sequence.length);
                del_symbols((&S1_set->sequences[i]->sequence), glob_period_end, S1_set->sequences[i]->sequence.length);
            }

            if (S_real_begin < glob_period_begin) {
                del_symbols((&S_set->sequences[i]->sequence), 0, glob_period_begin);
                del_symbols((&S1_set->sequences[i]->sequence), 0, glob_period_begin);
            } else {
                del_symbols((&S_set->sequences[i]->sequence), 0, S_real_begin);
                del_symbols((&S1_set->sequences[i]->sequence), 0, S_real_begin);
            }
        } else {
            S_set->sequences[i]->orig_begin = 0;
            S1_set->sequences[i]->orig_begin = 0;
            S_set->sequences[i]->orig_end = seq_set->sequences[i]->length;
            S1_set->sequences[i]->orig_end = (size_t)period_sym;
        }

        for (long j = S_set->sequences[i]->sequence.length - 1; j >= 0; --j) {
            if (S_set->sequences[i]->sequence.symbols[j] >= S_set->sequences[i]->sequence.count_symbols &&
                S1_set->sequences[i]->sequence.symbols[j] >= S1_set->sequences[i]->sequence.count_symbols
            ){
                del_symbol(&(S_set->sequences[i]->sequence), j);
                del_symbol(&(S1_set->sequences[i]->sequence), j);
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    double R = sqrt(calculate_R_sqr(seq_set->sequences[0]->count_symbols, (size_t)period_sym, config.R_sqr_mult_const));
    
    mahds_matrix* freq_matrix = mahds_matrix_init(seq_set->sequences[0]->count_symbols, (size_t)period_sym);
	create_frequency_matrix(S_set, S1_set, freq_matrix);
	mahds_matrix* weight = mahds_matrix_init(freq_matrix->lines, freq_matrix->columns);
	create_weight_matrix(freq_matrix, weight);
	mahds_matrix* normalized_weight = mahds_matrix_init(weight->lines, weight->columns);
	normalize_weight_matrix(weight, seq_set, R, config.K_d, (size_t)period_sym, normalized_weight);
    mahds_matrix_free(freq_matrix);
    mahds_matrix_free(weight);

    double F_max = 0.;
    char gapS = 0;
    char gapS1 = 0;
    for (size_t i = 0; i < S_set->seq_number; ++i) {
        for (size_t j = 0; j < S_set->sequences[i]->sequence.length; ++j) {
            if (S_set->sequences[i]->sequence.symbols[j] >= S_set->sequences[i]->sequence.count_symbols) {
                if (gapS) {
                    F_max -= config.e_price;
                } else {
                    F_max -= config.d_price;
                }
                gapS = 1;
                gapS1 = 0;
            } else if (S1_set->sequences[i]->sequence.symbols[j] >= S1_set->sequences[i]->sequence.count_symbols) {
                if (gapS1) {
                    F_max -= config.e_price;
                } else {
                    F_max -= config.d_price;
                }
                gapS = 0;
                gapS1 = 1;
            } else {
                F_max += mahds_matrix_get(normalized_weight, S_set->sequences[i]->sequence.symbols[j], S1_set->sequences[i]->sequence.symbols[j]);
                gapS = 0;
                gapS1 = 0;
            }
        }
    }
    
    
    for (size_t i = 0; i < seq_set->seq_number; ++i) {
        free(S_set->sequences[i]->sequence.symbols);
        free(S_set->sequences[i]);

        free(S1_set->sequences[i]->sequence.symbols);
        free(S1_set->sequences[i]);

        free(seq_set->sequences[i]->symbols);
        free(seq_set->sequences[i]);
    }
    free(S_set->sequences);
    free(S_set);
    free(S1_set->sequences);
    free(S1_set);
    free(seq_set->sequences);
    free(seq_set);

    for (size_t i = 0; i < seq_inf->seq_set->seq_number; ++i) {
        free(seq_inf->seq_set->sequences[i]->symbols);
        free(seq_inf->seq_set->sequences[i]);
        free(seq_inf->mapping_info->original_headers[i]);
        free(seq_inf->mapping_info->original_set[i]->symbols);
        free(seq_inf->mapping_info->original_set[i]);
    }
    free(seq_inf->seq_set->sequences);
    free(seq_inf->seq_set);
    free(seq_inf->mapping_info->original_set);
    free(seq_inf->mapping_info->original_headers);
    free(seq_inf->mapping_info);
    free(seq_inf);


    mahds_matrix_free(normalized_weight);
    free(config.matrix_path);
    
    return F_max;
}
