#!/usr/bin/python

import os
import subprocess
import datetime
import re
import sys
sys.path.append("..")
from common.common import (
    get_config_file,
    get_json_data,
)

re_bad_stdout = re.compile(r"(?:invalid)|(?:error)", flags=re.IGNORECASE)

__path = os.path.dirname(os.path.abspath(__file__))
__default_config_file = __path+r"/eet_config.json"


def process_all_fams(src_dirname:str, dst_dirname:str, client_script:str, cl_params:dict, log_file:str) -> None:
    if not os.path.exists(src_dirname) or not os.path.isdir(src_dirname):
        raise Exception(f'"{src_dirname}" is not a directory.')
    if not os.path.exists(dst_dirname) or not os.path.isdir(dst_dirname):
        os.mkdir(dst_dirname)

    for fname in os.listdir(src_dirname):
        if fname.find(".fasta") == -1 and fname.find(".tfa") == -1:
            continue
            
        full_fname = os.path.join(src_dirname, fname)
        output_fname = ""
        output_fname = fname.replace(".tfa",".txt")
        output_fname = output_fname.replace(".fasta",".txt")
        full_output_fname = os.path.join(dst_dirname, output_fname)

        if not cl_params['rewrite'] and os.path.exists(output_fname) and not os.path.isdir(output_fname):
            continue
        
        method = os.path.basename(client_script).replace(".py","")
        title = cl_params["method_params"][method]["out_prefix"] + output_fname.replace(".txt","")
        execution_params_str = f"python {client_script} " + cl_params["method_params"][method]["base_param_str"] + f" --outfile {title} " + full_fname

        executed = subprocess.Popen(
            execution_params_str.split(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding='utf-8'
        )
        return_code = executed.wait()
        stdout, stderr = executed.communicate()

        with open(log_file, 'a') as log:
            t = datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S").rstrip()
            msg = "NORMALLY EXECUTED"
            if re_bad_stdout.search(stdout):
                msg = "ERROR(?) OCCURRED"
            if return_code or stderr:
                msg = "ERROR OCCURRED"
            log.write(f">>> status: [{msg}]\ndatetime: [{t}]\nscript: [{execution_params_str.rstrip()}]:\nreturn code: [{return_code}]\nstdout: [{stdout.rstrip()}]\nstderr: [{stderr.rstrip()}]\n")

        cwd = os.getcwd()
        aln_file = ""
        other_generated_files = []
        re_aln_file = re.compile(f'^{title}'+r'\.aln.*?fasta\.fasta$')

        for cur_fname in os.listdir(cwd):
            if re_aln_file.search(cur_fname):
                aln_file = cur_fname
            else:
                if cur_fname.find(title) == 0:
                    other_generated_files.append(cur_fname)
        
        if aln_file:
            os.replace(os.path.join(cwd, aln_file), full_output_fname)
        # else:
            # raise Exception('Alignment file not found')

        for f in other_generated_files:
            os.remove(os.path.join(cwd, f))


def launch_alignments() -> None:
    config = get_json_data(get_config_file(__default_config_file))

    clients_dir = config["clients_dir"]
    if not os.path.exists(clients_dir) or not os.path.isdir(clients_dir):
        raise Exception(f'"{clients_dir}" is not a directory.')
    os.chdir(clients_dir)
    for sname in os.listdir(clients_dir):
        if sname.find(".py") == -1:
            continue
        fullsname = os.path.join(clients_dir, sname)

        if not os.path.exists(config["dst_dir"]) or not os.path.isdir(config["dst_dir"]):
            os.mkdir(config["dst_dir"])

        res_dir = os.path.join(config["dst_dir"], sname.replace(".py",""))
        if not os.path.exists(res_dir) or not os.path.isdir(res_dir):
            os.mkdir(res_dir)

        process_all_fams(config["src_dir"], res_dir, fullsname, config["params"], config["log_file"])



if __name__ == "__main__":
    launch_alignments()

# python clustalo.py --email dk0stenko@yandex.ru --stype protein --outfmt fa --order input --outfile TestAlign --quiet ~/Desktop/BaliBase/BAliBASE_R10/RV100/BBA0013.tfa

# python kalign.py --email dk0stenko@yandex.ru --stype protein --format fasta --outfile TestAlign --quiet ~/Desktop/BaliBase/BAliBASE_R10/RV100/BBA0013.tfa

# python mafft.py --email dk0stenko@yandex.ru --stype protein --format fasta --order input --outfile TestAlign --quiet ~/Desktop/BaliBase/BAliBASE_R10/RV100/BBA0013.tfa

# python muscle.py --email dk0stenko@yandex.ru --format fasta --outfile TestAlign --quiet ~/Desktop/BaliBase/BAliBASE_R10/RV100/BBA0013.tfa

# python prank.py --email dk0stenko@yandex.ru --outfile TestAlign --quiet ~/Desktop/BaliBase/BAliBASE_R10/RV100/BBA0013.tfa
# TestAlign.aln-fasta.fasta // also TestAlign.aln-1-fasta.fas TestAlign.aln-2-fasta.fas

# python tcoffee.py --email dk0stenko@yandex.ru --stype protein --format fasta_aln --order input --outfile TestAlign --quiet ~/Desktop/BaliBase/BAliBASE_R10/RV100/BBA0013.tfa