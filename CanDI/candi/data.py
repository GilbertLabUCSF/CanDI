# data.py loads data automatically and contains methods for loading data via user input
import os
import gc
import operator
import json
import configparser
from pathlib import Path
import pandas as pd
import numpy as np
import sys

class Data(object):
    """Class data is used for loading and caching data
    can be tuned to load specific datasets upon import by editing config.ini
    can call Data.load() to load any specific dataset
    """
    def __init__(self):

        self._file_path = Path(os.path.dirname(os.path.realpath(__file__))).parent.absolute()
        config_path = self._file_path / 'setup/data/config.ini'
        print(self._file_path)
        print(config_path)

        parser = configparser.ConfigParser() #parses config for data sources
        parser.read(config_path)

        self._parser = parser
        self._verify_install()
        self._init_sources()
        self._init_depmap_paths()
        self._init_index_tables()

    def _verify_install(self): #ensures data being loaded is present
        try:
            assert "depmap_urls" in self._parser.sections()
        except AssertionError:
            raise(RuntimeError, "CanDI has not been properly installed. Please run CanDI/install.py prior to import")

    def _init_sources(self):
        """this function creates paths
        to source directories of each data source."""
        sources = []

        for option in self._parser["data_paths"]:
            sources.append(option)
            setattr(self, "_" + option + "_path", self._file_path / self._parser["data_paths"][option])

        self.sources = sources

    def _init_depmap_paths(self):

        """This function establishes file paths
        for datasets retrieved from depmap."""

        self.depmap_files = self._parser["depmap_files"]

        for option in self._parser["depmap_files"]:
            try:
                new_path = self._depmap_path / self._parser.get("depmap_files", option)
                assert os.path.exists(new_path)
                setattr(self, option, new_path)
            except AssertionError:
                setattr(self, option, None)

    def _init_index_tables(self):

        for option in self._parser["autoload_info"]:
            try:
               new_path = self._file_path / self._parser.get("autoload_info", option)
               print(new_path)
               assert os.path.exists(new_path)
               setattr(self, option, self._handle_autoload(option, new_path))
            except AssertionError:
               raise(RuntimeError, "You are missing essential index table: {}. exiting...".format(new_path))


    @staticmethod
    def _handle_autoload(method, path):
        """This function loads datasets
        that are under the autoload section
        of the config.ini file.
        """
        if method == "genes":

            df = pd.read_csv(path,
                             memory_map=True,
                             low_memory=False,
                             dtype={"ENTREZ ID":str},
                             index_col = "Approved symbol")

        elif method == "cell_lines":

            df = pd.read_csv(path,
                             memory_map=True,
                             low_memory=False,
                             index_col="DepMap_ID")

        elif method == "locations":
            df = pd.read_csv(path,
                             memory_map=True,
                             low_memory=True)

        return df


    def load(self, key):
        """This function loads a dataset into memory as a pandas DataFrame.

        Note:
            If the nth row and mth column is equal to 0 the nth gene is not mutated in the mth cell line.
            If the nth row and mth column is equal to 1 the nth gene is mutated in the mth cell line.
        Args:
            key:
        Returns:
            [---]
        """
        if hasattr(self, key):

            new_path = getattr(self, key)
            try:
                index = self._parser.get("index", key)

            except configparser.NoOptionError:
                index = None

            except AttributeError:
                raise RuntimeError("CanDI is not compatible with python2. Please ensure you're using Python3.")

            df = pd.read_csv(new_path,
                             memory_map = True,
                             low_memory = False,
                             index_col = index)

            setattr(self, key, df)
            return getattr(self, key)

        else:
            raise KeyError("{0} cannot find file {1}".format(self, key))


    def unload(self, key):
        """This function removes a dataset from memory
        Args:
            key:
        """
        setattr(self, key, self.data_path + self._parser["files"][key])
        gc.collect()
