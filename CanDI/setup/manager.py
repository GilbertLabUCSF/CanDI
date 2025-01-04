"""
The manager module handles interations with the datasources
and the config file. It is used to setup of the config file upon installation.
All data downloading is done by Manager class and its subclasses.
"""

import os
import configparser
import json
import time
import requests
import numpy as np
import polars as pl
import pandas as pd
from time import sleep
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from . import dataverse


class Manager(object):
    def __init__(self, manager_path='auto', cfig_path='auto', verbose=False):
        """Initializes the Manager class

        Args:
            manager_path (str, optional): The path to the manager directory. This is where the data will be stored.
            cfig_path (str, optional): The path to the config file.
        """
        if manager_path == 'auto':
            manager_path = os.path.dirname(os.path.realpath(__file__))
        else:
            # make sure the path is a directory and exists or create it
            if not os.path.exists(manager_path):
                os.makedirs(manager_path)

        if cfig_path == 'auto':
            if not os.path.exists(manager_path + "/data/config.ini"):
                # copy config.draft.ini to config.ini
                os.system(f"cp {manager_path}/data/config.draft.ini {manager_path}/data/config.ini")
            cfig_path = manager_path + "/data/config.ini"

        elif type(cfig_path) == str and os.path.exists(cfig_path):
            cfig_path = cfig_path
        
        elif type(cfig_path) == str and not os.path.exists(cfig_path):
            raise FileNotFoundError(f"Config file not found at {cfig_path}")

        if verbose:
            print(f"Manager Path: {manager_path}")
            print(f"Config Path: {cfig_path}")

        parser = configparser.ConfigParser()
        parser.read(cfig_path)

        self.manager_path = manager_path
        self.cfig_path = Path(cfig_path)
        self.parser = parser

    @staticmethod    
    def write_config(cfig_path, parser):

        print("Writing config file")
        with open(cfig_path, "w") as f:
            parser.write(f)
            f.close()


class DataverseDepMap(Manager):
    def __init__(self, manager_path='auto', cfig_path='auto', verbose=False):
        super().__init__(manager_path, cfig_path, verbose)
        self.release = '21Q4' # default release uploded to CanDI dataverse
        self.download_source = 'dataverse, ' + dataverse.CANDI_DATAVERSE_DOI
    
    def download_reformatted_data(self):
        if not os.path.exists(self.manager_path + '/data/'):
            os.makedirs(self.manager_path + '/data/')

        if not os.path.exists(self.manager_path + '/data/depmap/'):
            os.makedirs(self.manager_path + '/data/depmap/')

        session = dataverse.DepMapDownloader()
        urls, file_names = session.download(
            self.manager_path + '/data/depmap/', 
            return_type= ["url", "name"]
        )

        depmap_urls = {
            file: url for url, file in zip(urls, file_names)
        }

        depmap_files = {}
        for file in file_names:
            f_key = file.split('.')[0]
            f_key = f_key.replace('CCLE_','')
            f_key = f_key.replace('CRISPR_','')
            depmap_files[f_key] = file 

        formatted = {
            f'{self.manager_path}/data/depmap/{file}': file for file in file_names 
            if 'readme' not in file.lower()
        }

        data_paths = {
            'depmap':   'data/depmap/',
            'genes':    'data/genes/',
            'corum':    'data/complexes/',
            'location': 'data/location/'
        }
        
        self.parser["depmap_urls"] = depmap_urls
        self.parser["depmap_files"] = depmap_files
        self.parser["formatted"] = formatted
        self.parser["data_paths"] = data_paths


class BroadDepMap(Manager):
    def __init__(self, manager_path='auto', cfig_path='auto', verbose=False):
        super().__init__(manager_path, cfig_path, verbose)
        self.download_source = 'Broad DepMap, https://depmap.org/'
        
    def get_depmap_info(self, release="latest"):

        depmap = self.parser["download_urls"]["depmap"]
        print("Getting download information from DepMap")
        response = requests.get(depmap)
        assert response.status_code == 200
        print("GET Successful")

        self.response = response.json()
        self.release = self.get_release(release)
        self.download_info, self.depmap_files = self.parse_release()
        self.parser["depmap_urls"] = self.download_info
        self.parser["depmap_files"] = self.depmap_files

    def parse_release(self):

        download_urls = {}
        depmap_files = {}
        for table in self.response["table"]:

            if self.release == table["releaseName"] and table["downloadUrl"]:

                download_urls[table["fileName"]] = table["downloadUrl"]
                depmap_files[self.format_filename(table["fileName"])] = table["fileName"]

        return download_urls, depmap_files

    def get_release(self, release):

        if release == "latest":
            release_info = [i for i in self.response["releaseData"] if i["isLatest"] is True][0]

        else:
            release_info = [i for i in self.response["releaseData"] if release in i["releaseName"]][0]

        self.parser["depmap_release"] = release_info

        return release_info["releaseName"]

    def format_filename(self, filename, release):

        # set candi_name to the filename without the extension
        candi_name = filename.split(".")[0]

        if release == "21Q4":
            if "CRISPR_" in candi_name:
                candi_name = candi_name[len("CRISPR_"):]
            elif "CCLE_" in candi_name:
                candi_name = candi_name[len("CCLE_"):]
            if 'v2' in candi_name:
                candi_name = candi_name[:-len("_v2")]
        else:
            #TODO: add more cases for different releases, e.g. 24Q4 new file formats
            pass

        return candi_name

    def depmap_download(self, name, filename=False):

        time.sleep(1)
        entry = self.manage_request(name, "depmap")
        self.fetch_url(entry)

    def fetch_url(self, entry):

        filename, path, url = entry

        r = requests.get(url, stream=True)
        length = int(r.headers["content-length"])
        chunk_size = 1024
        total = 0

        print("Downloading {}...".format(filename))
        if r.status_code == 200:

            with open(path, "w") as f:

               for chunk in r.iter_content(chunk_size=chunk_size):

                   total += chunk_size
                   chunk = chunk.decode("utf-8")
                   chunk = chunk.replace("\t", ",")
                   f.write(chunk)
            print("Downloading {} complete!".format(filename))
            downloads = self.parser["downloads"]

            downloads[filename] = str(path)

    def parallel_fetch(self, entries):
        print("Starting Pool")
        with ThreadPoolExecutor(max_workers=4) as executor:
            for i in entries:
                executor.submit(self.fetch_url, i)

    def download_defaults(self):

        default_sources = json.loads(self.parser.get("defaults","downloads"))
        to_download =  json.loads(self.parser.get("defaults", default_sources[0])) 

        entries = [self.manage_request(i, "depmap") for i in to_download]
        self.parallel_fetch(entries)

    def manage_request(self, name, path, filename=False):

        if filename:
            filename = name
        else:
            filename = self.parser['depmap_files'][name]

        url = self.parser['depmap_urls'][filename]

        base_path = Path(self.manager_path) / self.parser["data_paths"][path]

        try:
            assert os.path.exists(base_path)
        except AssertionError:
            os.mkdir(base_path)

        write_path = base_path / filename

        return (filename, write_path, url)

    def depmap_autoformat(self):

        try:
            downloaded = self.parser["downloads"] 
        except KeyError:
            raise(RuntimeError, "There are not data files to format. Please download data and try again or run install.py")

        for k,v in downloaded.items():
            try:
                if k not in self.parser["formatted"]:
                    print("Formatting {}".format(k))
                    df = pd.read_csv(v, low_memory=False, memory_map=True)
                    self.format_depmap_data(df, v)
            except KeyError:
                print("Formatting {}".format(k))
                df = pd.read_csv(v, low_memory=False, memory_map=True)
                self.format_depmap_data(df, v)

    def format_depmap_data(self, df, path, release):

        if release == "21Q4":
            if ("AAAS (8086)" in df.columns) or ("AAAS (ENSG00000094914)" in df.columns):

                df.rename(columns = lambda s: s.split(" ")[0], inplace=True)

                if "Unnamed:" in df.columns:
                    df.rename(columns={"Unnamed:":"DepMap_ID"}, inplace=True)

                df = df.set_index("DepMap_ID").T
                df.reset_index(inplace=True)
                df.rename(columns={"index":"gene"}, inplace=True)
                df.set_index("gene", inplace=True)
                df.to_csv(path)

            if "Protein_Change" in df.columns:

                try:
                    df.drop("Unnamed: 0", axis=1, inplace=True)
                    df.to_csv(path, index=False)
                except KeyError:
                    pass

            if "Hugo_Symbol" in df.columns:
                try:
                    df.rename(columns={"Hugo_Symbol": "gene"}, inplace=True)
                    df.to_csv(path, index=False)
                except KeyError:
                    pass

            if "LeftGene" in df.columns:
                for col in df.columns:
                    if "Gene" in col:
                        split_cols = df[col].str.split(" ", expand=True)
                        df[col] = split_cols[0]
                        df[col[:-4] + "EnsemblID"] = split_cols[1].str.replace("(", "").str.replace(")", "")

                df.to_csv(path, index=False)
        else:
            #TODO: add more cases for different releases, e.g. 24Q4 new file formats
            pass

        try:
            formatted = self.parser["formatted"]
        except KeyError:
            self.parser["formatted"] = {}
            formatted = self.parser["formatted"]

        formatted[path.split("/")[-1]] = path


class SangerDepMap(Manager):
    def __init__(self, cfig_path='auto'):
        super().__init__(cfig_path)

    def sanger_download():
        pass


class DataverseCoessentiality(Manager):
    def __init__(self, manager_path='auto', cfig_path='auto', verbose=False):
        super().__init__(manager_path, cfig_path, verbose)
        self.download_source = 'Dataverse'
        self.reference = 'https://github.com/kundajelab/coessentiality'
        self.verbose = verbose
    
    def download_raw_files(self):
        if not os.path.exists(self.manager_path + '/data/'):
            os.makedirs(self.manager_path + '/data/')

        if not os.path.exists(self.manager_path + '/data/coessentiality/'):
            os.makedirs(self.manager_path + '/data/coessentiality/')
        
        session = dataverse.CoessentialityDownloader()
        urls, file_names = session.download(
            self.manager_path + '/data/coessentiality/',
            return_type= ["url", "name"]
        )

        self.urls = urls
        self.file_names = file_names

    def _load_coessentiality_matrix(self):
        data_dir = f'{self.manager_path}/data/coessentiality'

        gene_names = pd.read_csv(
            f'{data_dir}/genes.txt',header=None,names=['gene_name']
        )['gene_name']
        
        GLS_sign = np.load(f'{data_dir}/GLS_sign.npy')
        GLS_p = np.load(f'{data_dir}/GLS_p.npy')
        
        self.matrix = pl.from_dataframe(
            pd.DataFrame((-1*np.log10(GLS_p)) * GLS_sign, columns = gene_names, index = gene_names).reset_index()
        )

    def _get_coessentiality_df(self, pvalue_threshold = 10**-3):
        df = self.matrix.melt('gene_name')
        df.columns = ['gene_1','gene_2','coessentiality']
        df = df.filter(~(pl.col('gene_1') == pl.col('gene_2')))
        df = df.filter(pl.col('coessentiality') > -np.log10(pvalue_threshold))
        
        self.df = df
        self.pvalue_threshold = pvalue_threshold

    def coessentiality_autoformat(self):
        
        coessentiality_matrix_path = f'{self.manager_path}/data/coessentiality/coessentiality_matrix.csv'
        coessentiality_df_path = f'{self.manager_path}/data/coessentiality/coessentiality_df.csv'

        # Check if the data has already been formatted or run the formatting
        if os.path.exists(coessentiality_matrix_path):
            if self.verbose: print("coessentiality_matrix.csv already exists")
        
        else:
            if self.verbose: print("Building Coessentiality Matrix ...", end=' ')
            self._load_coessentiality_matrix()
            self.matrix.to_pandas().to_csv(coessentiality_matrix_path)
            if self.verbose: print("Done!")
        
        if os.path.exists(coessentiality_df_path):
            if self.verbose: print("coessentiality_df.csv already exists")
        
        else:
            if self.verbose: print("Building Coessentiality DataFrame ...", end=' ')
            self._get_coessentiality_df()
            self.df.to_pandas().to_csv(coessentiality_df_path)
            if self.verbose: print("Done!")
        
        # Update the config file
        self.parser['data_paths'].update({
            'coessentiality': 'data/coessentiality/'
        })

        self.parser['formatted'].update({
            'coessentiality_matrix.csv': coessentiality_matrix_path,
            'coessentiality_df.csv': coessentiality_df_path
        })
        
        self.parser['depmap_files'].update({
            'coessentiality': coessentiality_df_path,
            'coessentiality_matrix': coessentiality_matrix_path,
        })
