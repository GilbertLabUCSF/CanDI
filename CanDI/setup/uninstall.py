import os
import shutil
import argparse
from .manager import Manager


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", help="Specify the database to uninstall", default="depmap")
    parser.add_argument("--directory", help="Specify the data parent directory", default='auto')
    args = parser.parse_args()

    if args.database == 'depmap':
        print("Uninstalling DepMap data")

        m = Manager()

        if args.directory == 'auto':
            shutil.rmtree(m.manager_path + "/data/depmap/")
        elif os.path.exists(args.directory):
            shutil.rmtree(m.manager_path + "/data/depmap/")
        else:
            raise ValueError("Invalid data directory")
    else:
        raise ValueError("Invalid database. Currently only 'depmap' is supported")
    
if __name__ == "__main__":
    main()