import random
import os
from functools import reduce
from common.common import (
    get_config_file,
    get_json_data
)

canonical_aminoacids = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
class RndSequencesGenerator:
    def __init__ (self, sequence_length:int=0, sequences_number:int=1, symbol_types:int=20, symbol_first:str='A', symbols_set:tuple=None):
        self.sequences_number = sequences_number
        self.sequence_length = sequence_length
        self.symbol_types = symbol_types
        self.symbol_first = symbol_first
        self.symbols_set = symbols_set
        random.seed()

    def generate_sequence(self):
        res = [0]*self.sequence_length
        if not self.symbols_set:
            left_boud = ord(self.symbol_first)
            right_boud = left_boud + self.symbol_types - 1
            for i in range(len(res)):
                res[i] = chr(random.randint(left_boud, right_boud))
        else:
            r_idx = len(self.symbols_set)-1
            for i in range(len(res)):
                res[i] = self.symbols_set[random.randint(0,r_idx)]
        return res

    def generate_seq_set(self, seq_length_list:list=None):
        if seq_length_list:
            res = [None]*len(seq_length_list)
            base_seq_length = self.sequence_length
            for i in range(len(seq_length_list)):
                self.sequence_length = seq_length_list[i]
                res[i] = self.generate_sequence()
            self.sequence_length = base_seq_length
            return res   
        else:
            return [self.generate_sequence() for _ in range(self.sequences_number)]


def output_fasta_file(seq_set:list, file_path:str):
    with open(file_path, 'w') as file:
        for i in range(len(seq_set)):
            file.write("> {}\n".format(i) + "".join(seq_set[i]) + "\n")

#####################################################################################


__path = os.path.dirname(os.path.abspath(__file__))
__default_config_file = __path+r"/rsg_config.json"
def launch() -> None:
    cfg = get_json_data(get_config_file(__default_config_file))
    rsg = RndSequencesGenerator(symbols_set = canonical_aminoacids)
    seq_length_list = reduce(lambda a, b: a + [int(b["length"])] * int(b["number"]), cfg["seqs_conf"], [])
    
    if not os.path.exists(cfg["dst_dir"]) or not os.path.isdir(cfg["dst_dir"]):
        raise Exception(f"Destination dir ({cfg['dst_dir']}) does not exist.")
    
    seq_conf_str = reduce(lambda a, b: a + f"l{int(b['length'])}n{int(b['number'])}", cfg["seqs_conf"], "")
    fulldirname = os.path.join(cfg["dst_dir"], f"sets_{seq_conf_str}")
    if not os.path.exists(fulldirname) or not os.path.isdir(fulldirname):
        os.mkdir(fulldirname)

    for i in range(int(cfg["sets_num"])):
        random.shuffle(seq_length_list)
        output_fasta_file(
            rsg.generate_seq_set(seq_length_list),
            os.path.join(fulldirname,f"set{i}_{seq_conf_str}.fasta")
        )


if __name__ == "__main__":
    launch()
