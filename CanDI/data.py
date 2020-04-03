#data.py loads data automatically and contains methods for loading data via user input
import os
import gc
import operator
import configparser
import pandas as pd
import numpy as np
import sys



class Data(object):
#Class data is used for loading and caching data
#can be tuned to load specific datasets upon import by editing config.ini
#can call Data.load() to load any specific dataset
    def __init__(self):

        file_path = os.path.dirname(os.path.abspath(__file__))
        config_path = file_path + '/data/config.ini'
        parser = configparser.ConfigParser()
        parser.read(config_path)
        self._parser = parser
        #Path used for loading data
        data_path = file_path + parser["paths"]["datasource"]
        self.data_path = data_path
        #Assign File names to variables
        files = parser['files']
        info = parser["info"]

        #set attributes to files paths that are accesibile by data
        self.transcription = data_path + files["transcription"]
        self.counts = data_path + files["counts"]
        self.pickles = data_path + files["pickles"]
        self.complexes = data_path + files["complexes"]
        self.copy_number = data_path + files["copy_number"]
        self.mutations = data_path + files["mutations"]
        self.depmap = data_path + files["depmap"]
        self.fusions = data_path + files["fusions"]
        self.translocations = data_path + files["translocations"]
        self.interactions = data_path + files["protein_interactions"]
        self.methylation_p = data_path + files["methylation_p"]
        self.methylation_1kb = data_path +files["methylation_1kb"]

        #Data will always load genes and cell line info because they're required to instantiate gene_query classes
        self.genes = pd.read_csv(data_path+info["genes"],
                                 memory_map=True,
                                 low_memory=False,
                                 dtype={"ENTREZ ID":str},
                                 index_col = "Approved symbol")
        self.cell_lines = pd.read_csv(data_path+info["cell_lines"],
                                      memory_map=True,
                                      low_memory=False,
                                      index_col="DepMap_ID")
        self.locations = pd.read_csv(data_path+files["locations"], memory_map=True, low_memory=True)

#        self._datasets = {"pickles":Pickles,
#                          "transcription":Transcription,
#                          "copy_number":CopyNumber,
#                          "depmap":"DepMap",
#                          "mutations":Mutations,
#                          "fusions":Fusions,
#                          "Translocations":Translocations,
#                          "interactions":Interactions,
#                          "complexes":Complexes,
#                          "methylation_p": "MethP",
#                          "methylation_1kb":"MethKb"}

        if bool(parser['behavior']['autoload']) == True:
            to_load = parser['behavior']['to_load']

            if type(to_load) == str:
                to_load = [to_load]

            for i in to_load:
                try:
                    setattr(self,
                            i, pd.read_csv(getattr(self, i),
                                                            memory_map=True,
                                                            low_memory=False,
                                                            index_col=self._parser['index'][i]))
                except KeyError:
                    setattr(self,
                            i,
                            (pd.read_csv(getattr(self, i),
                                                          memory_map=True,
                                                                          low_memory=False)))
                except AttributeError:
                    raise RuntimeError("This package is not compatible with python2. Make sure you're using a Python3 interpreter")

    def load(self, key):
        if hasattr(self, key):
            try:
                setattr(self,
                        key,
                        pd.read_csv(getattr(self, key),
                                    memory_map=True,
                                    low_memory=False,
                                    index_col=self._parser['index'][key]))
            except KeyError:
                setattr(self,
                        key,
                        pd.read_csv(getattr(self, key),
                                                       memory_map=True,
                                                       low_memory=False))
            except AttributeError:
                raise RuntimeError("This package is not compatible with python2. Make sure you're using a Python3 interpreter")

            return getattr(self, key)

        else:
            raise KeyError("{0} cannot find file {1}".format(self, key))

    def unload(self, key):

        setattr(self, key, self.data_path + self._parser["files"][key])
        gc.collect()
