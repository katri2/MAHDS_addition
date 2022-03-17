import os
import sys
import random
import copy
sys.path.append("..")
from rnd_sequences.rnd_seq_gen import RndSequencesGenerator, output_fasta_file
from common.common import (
    get_config_file,
    get_json_data
)

__path = os.path.dirname(os.path.abspath(__file__))
__default_config_file = __path+r"/sg_config.json"
canonical_aminoacids = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y") # canonical

class InheritedSequencesGenerator:
    def __init__ (
        self,
        symbols_set:tuple,
        pivot_seq:list,
        substitutions_cnt:int=0,
        insertions_cnt:int=0,
        deletions_cnt:int=0,
        substitutions_length:int=1,
        insertions_length:int=1,
        deletions_length:int=1
    ):
        self.pivot_seq = pivot_seq
        self.substitutions_cnt = substitutions_cnt
        self.insertions_cnt = insertions_cnt
        self.deletions_cnt = deletions_cnt
        self.substitutions_length = substitutions_length
        self.insertions_length = insertions_length
        self.deletions_length = deletions_length
        self.rnd_seq_gen = RndSequencesGenerator(symbols_set=symbols_set)

    def _do_rnd_insertion(self, res_seq:list) -> None:
        pos = random.randint(0, len(res_seq))
        self.rnd_seq_gen.sequence_length = self.insertions_length
        res_seq[pos:pos] = self.rnd_seq_gen.generate_sequence()

    def _do_rnd_deletion(self, res_seq:list) -> None:
        pos = random.randint(0, len(res_seq) - self.deletions_length)
        res_seq[pos:pos+self.deletions_length] = []

    def _do_rnd_substitution(self, res_seq:list) -> None:
        pos = random.randint(0, len(res_seq) - self.substitutions_length)
        self.rnd_seq_gen.sequence_length = self.substitutions_length
        res_seq[pos:pos+self.substitutions_length] = self.rnd_seq_gen.generate_sequence()

    def generate(self, order:tuple=None): #if None, then random. Other vars: ("s","d","i"), ("d", "i", "s"), ...
        res = copy.deepcopy(self.pivot_seq)
        action_seq = []
        if order:
            n = 0
            for action in order:
                if action == "i":
                    n = self.insertions_cnt
                elif action == "d":
                    n = self.deletions_cnt
                elif action == "s":
                    n = self.substitutions_cnt
                else:
                    raise Exception(f"In function 'generate(self, order=None)' unrecognized action: {action} in 'order'")
                action_seq += ([action] * n)
        else:
            action_seq = (["i"] * self.insertions_cnt) + (["d"] * self.deletions_cnt) + (["s"] * self.substitutions_cnt)
            random.shuffle(action_seq)

        for action in action_seq:
            if action == "i":
                self._do_rnd_insertion(res)
            elif action == "d":
                self._do_rnd_deletion(res)
            elif action == "s":
                self._do_rnd_substitution(res)
            else:
                raise Exception(f"In function 'generate(self, order=None)' unrecognized action: {action}")
        
        return res


def generate_from_config() -> None:
    config = get_json_data(get_config_file(__default_config_file))
    sg = RndSequencesGenerator(
        sequence_length = config["pivot_seq_length"],
        symbols_set = canonical_aminoacids
    )
    pivot_seq = sg.generate_sequence()
    
    if not os.path.exists(config["dst_dir"]) or not os.path.isdir(config["dst_dir"]):
        raise Exception(f"Destination dir ({config['dst_dir']}) does not exist.")
    
    isg = InheritedSequencesGenerator(canonical_aminoacids, pivot_seq)
    for params in config["params"]:
        isg.insertions_cnt = int(params["ic"])
        isg.insertions_length = int(params["il"])
        isg.deletions_cnt = int(params["dc"])
        isg.deletions_length = int(params["dl"])
        isg.substitutions_cnt = int(params["scm"] * config["pivot_seq_length"])
        isg.substitutions_length = int(params["sl"])
        
        fname = f"artif_ic{isg.insertions_cnt}il{isg.insertions_length}dc{isg.deletions_cnt}dl{isg.deletions_length}scm{params['scm']}sl{isg.substitutions_length}.fasta"
        fpath = os.path.join(config["dst_dir"], fname)
        
        res_set = [None] * config["heirs_num"]
        for i in range(config["heirs_num"]):
            res_set[i] = isg.generate()
        
        output_fasta_file(res_set, fpath)


if __name__ == "__main__":
    generate_from_config()
