import sys
import json


def get_config_file(default:str=None) -> str:
    if len(sys.argv) > 2:
        raise Exception(f'expected 0 or 1 command line argument, got: {len(sys.argv)-1}.')
    elif len(sys.argv) == 2:
        return sys.argv[1]
    else:
        return default


def get_json_data(fname:str) -> dict:
    with open(fname, "r") as file:
        return json.load(file)
