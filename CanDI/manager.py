import configparser
import requests
import time
import csv
import os


class Manager(object):


    def __init__(self):

        cfig_path = os.path.dirname(os.path.realpath(__file__)) + "/data/config.ini"
        parser = configparser.ConfigParser()
        parser.read(cfig_path)

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

    def from_depmap(self, filename):

        time.sleep(1)
        url = self.download_info[filename]

        print("Getting {}".format(filename))
        response = requests.get(url)
        text = response.content.decode("utf-8")

        with open(filename, "w", encoding="utf-8") as f:

            print("Writting {} to file".format(filename))
            for line in text.split("\n"):
                if line:
                    f.write(line)

            f.close()

    def write_config(self):

        print("Writing config file")
        with open(self.cfig_path, "w") as f:
            self.parser.write(f)
            f.close()

if __name__ == "__main__":

    m = Manager()
    m.get_depmap_info()
    m.write_config()
