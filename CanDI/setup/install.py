import argparse
from .manager import Manager

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", help="Specify the download source", default="dataverse")
    parser.add_argument("--data_dir", help="Specify the data directory", default=None)
    args = parser.parse_args()

    if args.source == 'dataverse':
        print("Downloading data from Dataverse")
        m = Manager(download_source=args.source, data_dir=args.data_dir)
        m.download_reformatted_data()
        m.write_config(m.cfig_path, m.parser)
    
    elif args.source == 'depmap':        
        print("Downloading data from DepMap")
        m = Manager(download_source=args.source, data_dir=args.data_dir)
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