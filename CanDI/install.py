from manager import Manager

if __name__ == "__main__":
    m = Manager()
    m.get_depmap_info()
    m.write_config(m.cfig_path, m.parser)
    m.download_defaults()
    m.write_config(m.cfig_path, m.parser)
    m.depmap_autoformat()
    m.write_config(m.cfig_path, m.parser)
