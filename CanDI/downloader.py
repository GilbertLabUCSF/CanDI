import configparser
import requests
import time
import csv

class Downloader(object):


    def __init__(self, release="latest"):

        parser = configparser.ConfigParser()
        parser.read("config.ini")

        url = parser["api"]["url"]

        print("Getting download information from DepMap")
        response = requests.get(url)
        assert response.status_code == 200
        print("GET Successful")

        self.response = response.json()
        self.parser = parser
        self.release = self.get_release(release)
        self.download_info = self.parse_release()


    def parse_release(self):

        download_urls = {}
        for table in self.response["table"]:

            if self.release == table["releaseName"]:

                download_urls[table["fileName"]] = table["downloadUrl"]

        return download_urls


    def get_release(self, release):

        if release == "latest":
            release_info = [i for i in self.response["releaseData"] if i["isLatest"] is True][0]

        else:
            release_info = [i for i in self.response["releaseData"] if release in i["releaseName"]][0]

        return release_info["releaseName"]

    def download(self, filename):

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

if __name__ == "__main__":

    d = Downloader()
    print(d.download_info)