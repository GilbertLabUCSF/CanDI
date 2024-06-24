import shutil
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", help="Specify the database to uninstall", default="depmap")
    parser.add_argument("--data_dir", help="Specify the data directory", default='auto')
    args = parser.parse_args()

    if args.database == 'depmap':
        print("Uninstalling DepMap data")
        shutil.rmtree(args.data_dir + "/data/depmap/")

    else:
        raise ValueError("Invalid database. Currently only 'depmap' is supported")
    
if __name__ == "__main__":
    main()