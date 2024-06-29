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
import subprocess


class Data(object):
    """Class data is used for loading and caching data
    can be tuned to load specific datasets upon import by editing config.ini
    can call Data.load() to load any specific dataset
    """
    def __init__(self, config_path='auto', verbose=False):

        if config_path == 'auto':
            self._file_path = Path(os.path.dirname(os.path.realpath(__file__))).parent.absolute() / 'setup'
            if os.path.exists(self._file_path / 'data/config.ini'):
                config_path = self._file_path / 'data/config.ini'
            else:
                config_path = self._file_path / 'data/config.draft.ini'
        
        elif os.path.exists(config_path) == False:
            raise FileNotFoundError("Config file not found at {}".format(config_path))
        elif os.path.exists(config_path) == True:
            if verbose: print("Using config file at {}".format(config_path))

        parser = configparser.ConfigParser() #parses config for data sources
        parser.read(config_path)

        self._parser = parser
        self._verify_install()
        self._init_sources()
        # self._init_depmap_paths()
        self._init_index_tables()

    def _verify_install(self): #ensures data being loaded is present
        #TODO: add more checks for different data sources
        try:
            assert "depmap_urls" in self._parser.sections()
        except AssertionError:
            print("CanDI has not been properly installed. ...")
            subprocess.run('candi-install', shell=True)

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
                assert os.path.exists(new_path)
                setattr(self, option, self._handle_autoload(option, new_path))
            except AssertionError:
                raise RuntimeError("You are missing essential index table: {}. exiting...".format(new_path))


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
                             sep='\t',
                             index_col="DepMap_ID")

        elif method == "locations":
            df = pd.read_csv(path,
                             memory_map=True,
                             low_memory=True)

        return df


    def load(self, key):
        """This function loads a dataset into memory as a pandas DataFrame.
        
        Args:
            key: str
               name of dataset to load into memory
        Returns:
            pandas.core.frame.DataFrame
                Pandas DataFrame is returned and saved as an attribute within the data object of the same name as key
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
            key: str
                name of the dataset to remove from memory
        Returns:
            None
                Data object of with same name as key is removed from memory and attribute is returned to dataset file path.
        """
        
        try:
            assert isinstance(getattr(self, key), pd.core.frame.DataFrame)
        except AssertionError:
            raise RuntimeError("{} is not currently loaded into memory".format(key))
            
        
        new_path = self._depmap_path / self._parser.get("depmap_files", key)
        assert os.path.exists(new_path)
        setattr(self, key, new_path)
        
