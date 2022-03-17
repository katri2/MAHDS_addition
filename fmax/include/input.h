#pragma once

#include "sequence.h"
#include "config.h"
#include "kseq.h"
#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>  

KSEQ_INIT(gzFile, gzread)


typedef enum {
    MAPPING_ARTIFICIAL,
    MAPPING_DNA,
    MAPPING_PROTEIN
} input_mapping;

typedef struct {
    size_t                   mapped_count_symbols;// how many symbol types present in fact (may be equal or more then mahds_symbolic_sequence->count_symbols)
    size_t                   original_count_symbols;
    char                     output_gap_symbol;// (in output seq set)
    char**                   original_headers;
    mahds_symbolic_sequence** original_set;
} mapping_info;

typedef struct {
    mahds_seq_set* seq_set;
    mapping_info* mapping_info;
} sequences_info;

sequences_info* fill_DNA_seq_info (kseq_t* original_set, char with_gaps); // with_gaps:bool
sequences_info* fill_protein_seq_info (kseq_t* original_set, char with_gaps); // with_gaps:bool
sequences_info* fill_artificial_seq_info (kseq_t* original_set, char with_gaps); // with_gaps:bool
sequences_info* parse_seq_info(char* file_path, input_mapping in_mapping, char with_gaps); // with_gaps:bool
void get_data_from_file(
    input_mapping         mapping_type,
    char*                 input_file_path,
    char*                 sys_cfg_path,
    char*                 align_cfg_path,
    mahds_config** OUT      config,
    sequences_info** OUT  seq_inf
);
