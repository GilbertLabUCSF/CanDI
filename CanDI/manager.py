import os
import sys
import configparser
import json
import time
from time import sleep
from pathlib import Path
import contextlib
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import requests

class Manager(object):
    """The Manager class handles interations with the datasources
    and the config file. It is used to setup of the config file upon installation.
    All data downloading is done by Manager
    """
    def __init__(self):

        manager_path = os.path.dirname(os.path.realpath(__file__))
        cfig_path = manager_path + "/data/config.ini"
        parser = configparser.ConfigParser()
        parser.read(cfig_path)

        self.manager_path = manager_path
        self.cfig_path = Path(cfig_path)
        self.parser = parser


    def sanger_download():
        pass

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

    def format_filename(self, filename):

        candi_name = filename.split(".")[0]

        if "Achilles_" in candi_name:
            candi_name = candi_name[len("Achilles_"):]
        elif "CCLE_" in candi_name:
            candi_name = candi_name[len("CCLE_"):]
        if 'v2' in candi_name:
            candi_name = candi_name[:-len("_v2")]

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

    def format_depmap_data(self, df, path):

        if "AAAS (8086)" in df.columns:

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

        try:
            formatted = self.parser["formatted"]
        except KeyError:
            self.parser["formatted"] = {}
            formatted = self.parser["formatted"]

        formatted[path.split("/")[-1]] = path


    @staticmethod
    def write_config(cfig_path, parser):

        print("Writing config file")
        with open(cfig_path, "w") as f:
            parser.write(f)
            f.close()

if __name__ == "__main__":
    m = Manager()
    #m.depmap_download("fusions")
    m.depmap_autoformat()
    m.write_config(m.cfig_path, m.parser)
