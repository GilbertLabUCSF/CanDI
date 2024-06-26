import os
import configparser
import json
from .manager import Manager


def write_cfig(cfig_path, parser):

    write_file = Manager.write_config
    write_file(cfig_path, parser)


def main(cfig_path='auto'):
    """
    This function will reset the config file to only contain the default sections.
    This is useful if you want to reset the config file to its original state.
    """
    if cfig_path == 'auto':
        cfig_path = os.path.dirname(os.path.realpath(__file__)) + "/data/config.ini"
    elif os.path.exists(cfig_path) == False:
        raise FileNotFoundError("Config file not found")
    elif os.path.exists(cfig_path) == True:
        print("Using config file at: " + cfig_path)
    
    parser = configparser.ConfigParser()
    parser.read(cfig_path)

    default_sections = json.loads(parser.get("defaults","sectionList"))

    print("Clearling non-default config sections")
    for sec in parser.sections():

        if sec not in default_sections:
            parser.remove_section(sec)

    write_cfig(cfig_path, parser)


if __name__ == "__main__":
    main()
