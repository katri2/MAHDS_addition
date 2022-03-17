from Bio import SeqIO  #sudo pip install Bio
import os
import sys


def get_fpath() -> str:
    if len(sys.argv) != 2:
        raise Exception(f'{sys.argv[0]} takes 1 CL parameter (path to file)')
    fpath = sys.argv[1]
    if not os.path.exists(fpath) or os.path.isdir(fpath):
        raise Exception(f'"{fpath}" must be a file with alignment.')
    return fpath

def count_gaps(seq, gap_sym:str="-") -> dict:
    """
    Returns gap openings count and overall gap count
    keys: openings, gaps
    """
    openings_and_gaps = {"openings": 0, "gaps": 0}
    prev_was_gap = False
    for sym in seq:
        if sym == gap_sym:
            openings_and_gaps["gaps"] += 1
            if not prev_was_gap:
                openings_and_gaps["openings"] += 1
                prev_was_gap = True
        else:
            prev_was_gap = False
    return openings_and_gaps

def count_gaps_fam(fasta_sequences, gap_sym:str="-") -> dict:
    """
    Returns gap openings count and overall gap count
    keys: openings, gaps
    """
    openings_and_gaps = {"openings": 0, "gaps": 0}
    for record in fasta_sequences:
        seq_res = count_gaps(record.seq)
        openings_and_gaps["openings"] += seq_res["openings"]
        openings_and_gaps["gaps"] += seq_res["gaps"]
    return openings_and_gaps
    
def launch() -> None:
    fpath = get_fpath()
    fasta_sequences = SeqIO.parse(fpath,'fasta')
    print(count_gaps_fam(fasta_sequences))


if __name__ == "__main__":
    launch()
