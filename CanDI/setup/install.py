import argparse
from .manager import DataverseDepMap, BroadDepMap


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", help="Specify the download source", default="dataverse")
    parser.add_argument("--data_dir", help="Specify the data directory", default='auto')
    args = parser.parse_args()

    if args.source == 'dataverse':
        print("Downloading data from Dataverse")
        m = DataverseDepMap(manager_path=args.data_dir, verbose=True)
        m.download_reformatted_data()
        m.write_config(m.cfig_path, m.parser)
    
    elif args.source == 'depmap':        
        print("Downloading data from DepMap")
        m = BroadDepMap(manager_path=args.data_dir, verbose=True)
        m.get_depmap_info()
        m.write_config(m.cfig_path, m.parser)
        m.download_defaults()
        m.write_config(m.cfig_path, m.parser)
        m.depmap_autoformat()
        m.write_config(m.cfig_path, m.parser)

    else:
        raise ValueError("Invalid source. Please specify either 'dataverse' or 'depmap'")
    
if __name__ == "__main__":
    main()