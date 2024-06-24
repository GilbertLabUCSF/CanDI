import os
import sys
import shutil
import argparse
from .manager import Manager


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", help="Specify the database to uninstall", default="depmap")
    parser.add_argument("--directory", help="Specify the data parent directory", default='auto')
    args = parser.parse_args()

    if args.database == 'depmap':
        print("Uninstalling CanDI: removing DepMap data")

        m = Manager()

        if args.directory == 'auto':
            depmap_path = m.manager_path + "/data/depmap/"
        elif os.path.exists(args.directory):
            depmap_path = args.directory + "/data/depmap/"
        else:
            sys.exit("Exit: Invalid directory path!")

        if not os.path.exists(depmap_path):
            sys.exit("Exit: Directory does not contain DepMap data")
        else:
            os.listdir(depmap_path)
            shutil.rmtree(depmap_path)
    else:
        raise ValueError("Invalid database. Currently only 'depmap' is supported")
    
if __name__ == "__main__":
    main()