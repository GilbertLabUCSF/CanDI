import numpy as np
from pathlib import Path
from CanDI.candi import data

class Grabber:
    """"Grabber class handles all bulk data retrival from the CanDI Classes.
    Behavior for data retrival is defined by the type of data being accessed
    as well as the class doing the retrival. Datasets are not all formatted the same way.
    """
    def __init__(self, grabber_type, key, axis):

        self.key = key
        self.gtype = getattr(self, grabber_type) #dict of with data as key and retrieval function as value
        self._isin_col = self.isin_dict[axis] #isin is a special type of subseting fuction
        self.axis = axis #different classes are indexed on different axes

    def __call__(self, item):
        """Core Grabber function, uses gtype dict to gather data.
        """
        if item not in dir(data):
            raise AttributeError("data has no attribute {}".format(item))

        dataset = getattr(data, item)
        if isinstance(dataset, Path):

            to_load = input("{} has not been loaded. Do you want to load, y/n?> ".format(item))
            if to_load == ("y" or "Y" or "Yes"):
                dataset = data.load(item)
                print("Load Complete")
            else:
                return

        return self.gtype[item](dataset)

    # """The following functions are the methods used for data retrival.
    # All datasets are loaded as pandas dataframes. These functions apply
    # standard pandas subsetting and indexing opperations.
    # """

    def get_one(self, dataset): #Get one element from user defined dataset

        cases = {0: lambda x,y: x.loc[y],
                 1: lambda x,y: x[y]}
        try:
            series = cases[self.axis](dataset, self.key)
        except KeyError:
            return

        return series

    def get_several(self, dataset): #Get several elements from user defined dataset

        getter = lambda x,y: x.reindex(y, axis=self.axis)
        values = getter(dataset, self.key).dropna(how="all", axis=self.axis)

        if values.empty:
            return
        else:
            return values

    def merge_two(self, dataset): #Merge two columns within one datasets

        getter = lambda x,y,z: x.loc[x[z].isin(y)]
        key = self.key
        try:
            assert type(key) is list
        except AssertionError:
            key = [key]

        left = getter(dataset, key, "LeftGene")
        right = getter(dataset, key, "RightGene")
        new_item = left.append(right).drop_duplicates()

        if new_item.empty:
            return
        else:
            return new_item

    def isin(self, dataset): #subset dataset if x is in a list of values
        key = self.key
        try:
            assert type(key) is list
        except AssertionError:
            key = [key]

        getter = lambda x, y, z: x.loc[x[z].isin(y)]
        item = getter(dataset, key, self._isin_col)
        if item.empty:
            return
        else:
            return item

    @property
    def isin_dict(self):

        return {0: "gene",
                1: "DepMap_ID"}
 
    # """The following Grabber properties are dictionaries with a dataset name as the key
    # and the retrieval function as the value. Each CanDI class has different op
    # """

    @property
    def gene(self): #defines data retrieval functions for gene class

        return {"gene_dependency": self.get_one,
                "gene_effect": self.get_one,
                "expression": self.get_one,
                "rnaseq_reads": self.get_one,
                "gene_cn": self.get_one,
                "locations": self.isin,
                "mutations": self.isin,
                "fusions": self.merge_two}
                #"complexes": self.get_one,
                #"translocations": self.merge_two,
                #"interactions": self.merge_two}
    @property
    def line(self): #defines data retrieval 

        return {"gene_dependency": self.get_one,
                "gene_effect": self.get_one,
                "expression": self.get_one,
                "rnaseq_reads": self.get_one,
                "gene_cn": self.get_one,
                "mutations": self.isin,
                "fusions": self.merge_two}
                #"translocations": self.isin}

    @property
    def canc(self):

        return {"gene_dependency": self.get_several,
                "gene_effect": self.get_several,
                "expression": self.get_several,
                "gene_cn": self.get_several,
                "mutations": self.isin,
                "fusions": self.isin,
                "rnaseq_reads": self.get_several}
                #"translocations": self.isin,

    @property
    def org(self):

        return {"gene_dependency": self.get_several,
                "gene_effect": self.get_several,
                "expression": self.get_several,
                "gene_cn": self.get_several,
                "rnaseq_reads": self.get_several,
                "mutations": self.isin,
                "fusions": self.merge_two}
                #"complexes": self.get_several,
                #"translocations": self.merge_two,
                #"interactions": self.merge_two}
