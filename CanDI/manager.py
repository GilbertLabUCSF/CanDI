import configparser
import requests
import time
import io
import csv
import os
import pandas as pd

class Manager(object):

    """
    The Manager class handles interations with the datasources
    and the config file. It is used to setup of the config file upon installation.
    All data downloading is done by Manager
    """

    def __init__(self):

        manager_path = os.path.dirname(os.path.realpath(__file__))
        cfig_path = manager_path + "/data/config.ini"
        parser = configparser.ConfigParser()
        parser.read(cfig_path)

        self.manager_path = manager_path
        self.cfig_path = cfig_path
        self.parser = parser

    def get_depmap_info(self, release="latest"):

        depmap = self.parser["download_apis"]["depmap"]
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
        elif 'v2' in candi_name:
            candi_name = candi_name[:-len("_v2")]

        return candi_name

    def depmap_download(self, filename):

        time.sleep(1)
        url = self.parser['depmap_urls'][filename]

        print("Downloading {}".format(filename))
        response = requests.get(url)
        content = response.content.decode('utf-8')

        print("Formatting {}".format(filename))
        df = pd.read_csv(io.StringIO(content))
        formatted = self.depmap_autoformat(df)

        self.write_df(filename, "depmap",formatted)

    def depmap_autoformat(self, df):

        if "AAAS (8086)" in df.columns:

            df.rename(columns = lambda s: s.split(" ")[0], inplace=True)
            df = df.set_index("Unnamed:").T
            df.reset_index(inplace=True)
            df.rename(columns={"index":"Gene"}, inplace=True)

        if "Protein_Change" in df.columns:
            df.drop("Unnamed: 0", axis=1, inplace=True)

        if "Hugo_Symbol" in df.columns:
            df.rename(columns={"Hugo_Symbol": "Gene"}, inplace=True)

        return df

    def write_df(self, filename, path, df):

        path = self.manager_path + self.parser["data_paths"][path]

        try:
            assert os.path.exists(path)
        except AssertionError:
            os.mkdir(path)

        print("Writting {0}".format(path+filename))
        df.to_csv(path+filename, index=False, sep=",")


        """
        with open(path+filename, "w", encoding="utf-8", newline='') as f:

            print("Writting {0}".format(path+filename))
            writer = csv.writer(f)
            writer.writerows(text)
            f.close()
        """

    @staticmethod
    def write_config(cfig_path, parser):

        print("Writing config file")
        with open(cfig_path, "w") as f:
            parser.write(f)
            f.close()

if __name__ == "__main__":

    m = Manager()
    m.depmap_download("ccle_mutations.csv")

