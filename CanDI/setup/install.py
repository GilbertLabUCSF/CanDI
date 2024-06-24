import argparse
from . import manager


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", help="Specify the database to download", default="depmap")
    parser.add_argument("--source", help="Specify the download source", default="dataverse")
    parser.add_argument("--directory", help="Specify the parent data directory", default='auto')
    args = parser.parse_args()

    if args.database == 'depmap':
        if args.source == 'dataverse':
            print("Downloading data from Dataverse")
            m = manager.DataverseDepMap(manager_path=args.directory, verbose=True)
            m.download_reformatted_data()
            m.write_config(m.cfig_path, m.parser)
        
        elif args.source == 'depmap':        
            print("Downloading data from DepMap")
            m = manager.BroadDepMap(manager_path=args.directory, verbose=True)
            m.get_depmap_info()
            m.write_config(m.cfig_path, m.parser)
            m.download_defaults()
            m.write_config(m.cfig_path, m.parser)
            m.depmap_autoformat()
            m.write_config(m.cfig_path, m.parser)

        else:
            raise ValueError("Invalid source. Please specify either 'dataverse' or 'depmap'")
    
    if args.database == 'coessentiality':
        if args.source == 'dataverse':
            print("Downloading data from Dataverse")
            m = manager.DataverseCoessentiality(manager_path=args.directory, verbose=True)
            m.download_raw_files()
            m.coessentiality_autoformat()
            m.write_config(m.cfig_path, m.parser)
        
        else:
            raise ValueError("Invalid source. Coessentiality data is only available on `dataverse`!")


if __name__ == "__main__":
    main()