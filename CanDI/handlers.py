import operator
import pandas as pd
import numpy as np
from pathlib import Path
import collections
import six
from CanDI import data


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
        """Core Grabber function, uses gtype dict to gather data."""
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

    """The following functions are the methods used for data retrival.
    All datasets are loaded as pandas dataframes. These functions apply
    standard pandas subsetting and indexing opperations.
    """

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
 
    """The following Grabber properties are dictionaries with a dataset name as the key
    and the retrieval function as the value. Each CanDI class has different op
    """

    @property
    def gene(self): #defines data retrieval functions for gene class

        return {#"pickles": self.get_one,
                "gene_effect": self.get_one,
                "expression": self.get_one,
                #"counts": self.get_one,
                "gene_cn": self.get_one,
                #"complexes": self.get_one,
                "locations": self.isin,
                "mutations": self.isin,
                "fusions": self.merge_two
                }
                #"translocations": self.merge_two,
                #"interactions": self.merge_two}
    @property
    def line(self): #defines data retrieval 

        return {#"pickles": self.get_one,
                "gene_effect": self.get_one,
                "expression": self.get_one,
                #"counts": self.get_one,
                "gene_cn": self.get_one,
                "mutations": self.isin,
                "fusions": self.merge_two
                #"translocations": self.isin}
            }


    @property
    def canc(self):

        return {#"pickles": self.get_several,
                "gene_effect": self.get_several,
                "expression": self.get_several,
                "gene_cn": self.get_several,
                "mutations": self.isin,
                "fusions": self.isin}
                #"translocations": self.isin,
                #"counts": self.get_several}

    @property
    def org(self):

        return {#"pickles": self.get_several,
                "gene_effect": self.get_several,
                "expression": self.get_several,
                "gene_cn": self.get_several,
                #"complexes": self.get_several,
                #"counts": self.get_several,
                "mutations": self.isin,
                "fusions": self.merge_two}
                #"translocations": self.merge_two,
                #"interactions": self.merge_two}

###################################################################################################

class BinaryFilter:
    """BinaryFilter class filters datasets based on a specific threshold.
    It's often useful to filter essentiality, expression, copy number etc.
    on specific thresholds. This class automates that behavior. BinaryFilter
    has different methods for handling different datatypes.
    """
    def __init__(self, margin, handler):

        self.margin = margin #number on which to filter
        self._default = self._handlers[handler]


    def __call__(self, vals, style, caller, threshold, return_lines=False):
        """Core function of binary filter. Behavior is defined during instantiation
        and is applied here.
        """

        self._eval_args(vals, style, threshold, return_lines)

        handler = self._handlers.get(type(vals), self._default)
        return handler(vals, style, caller, threshold, return_lines)


    def _float_handler(self, values, style, caller, *args):
        """This function handles filtering of numpy float objects.
        """

        if style == 'values':
            return values
        else:
            behaviors = {"over": operator.ge,
                         "under": operator.lt}

            return behaviors[caller](values, self.margin)


    def _series_handler(self, values, style, caller, *args):
        """This function handles filtering pandas series objects.
        """

        behaviors = {"over": values.ge,
                     "under": values.lt}

        evaluated = values[behaviors.get(caller)(self.margin)]

        if style == "values":
            return evaluated
        else:
            return list(evaluated.index)


    def _frame_handler(self, values, style, caller, threshold, return_lines):
        """This function handles filtering entire pandas dataframes.
        """

        behaviors = {"over": lambda x: x.gt(self.margin), 
                     "under": lambda x: x.lt(self.margin)}

        evaluated = values[behaviors.get(caller)(values)].dropna(thresh= int(threshold * values.shape[1]))

        if threshold != 1.0:
            evaluated = values.loc[evaluated.index]

        if style == "values":
            return evaluated

        elif return_lines is False:
            return list(evaluated.index)

        else:
            formatted = evaluated.to_dict("index").items()
            eva_dict =  {k:list(v.keys()) for k, v in formatted}
            return eva_dict


    @staticmethod
    def _eval_args(vals=None, style="bool", threshold=1.0, return_lines=False):
        """Function used to make sure arguments are correct.
        """

        assert style in ["bool", "values"],"style must be 'bool' or 'values'"
        assert 0.0 < threshold <= 1.0, "threshold is invalid, must be between 0 and 1"
        assert vals is not None

        try:
            assert ~vals.empty
        except AttributeError:
            pass
        return


    @property
    def _handlers(self): #Defines which handler to use based on the data type

        return {np.float64: self._float_handler,
                pd.Series: self._series_handler,
                pd.DataFrame: self._frame_handler}

###################################################################################################

class MutationHandler(object):
    """Handler filtering and querying of mutation data.
    Has methods for all CanDI objects that allow for more specific filtering.
    Includes methods for querying translocations and fusions.
    MutationHandler behavior is instantiated with instantiation of core CanDI objects.
    """

    def __init__(self, version):
        self.by = {"gene":"DepMap_ID",
                   "line":"gene",
                   "canc":"DepMap_ID",
                   "org": "DepMap_ID"}
        self.version = version


    def __call__(self, mut_dat, output, variant, item, translocations, fusions, all_except):
        """Core function of mutation handler.
        Behavior is defined on instantiation and applied in this function.
        """

        if variant and item:
            mut_dat = self._get_variant(mut_dat, variant, item, all_except=all_except) #get specific variant
        else:
            try:
                mut_dat = self._get_variant(mut_dat, "Variant_Classification", "Silent", all_except=True)
            except AssertionError:
                pass

        cases = {"gene": self._single_entity_mutated, #dict of functions to used based on version of handler being used
                "line": self._single_entity_mutated,
                "canc": self._multiple_entity_mutated,
                "org": self._multiple_entity_mutated}

        return cases[self.version](mut_dat, output, variant, item, translocations, fusions, all_except) #get mutations


    @staticmethod
    def _get_variant(mut_dat, variant, item, all_except=False):
        """Special case of mutation handler.
        Applied when user wants the specific variant
        of a specific mutation.
        """

        assert item in mut_dat[variant].unique(), "{0} not found, options are: {1}".format(item, mut_dat[variant].unique())


        if isinstance(item, collections.Iterable) and not isinstance(item, six.string_types):
            method = lambda x,y: mut_dat.loc[mut_dat[x].isin(y)]

        else:
            cases = {True: lambda x, y: mut_dat.loc[~(mut_dat[x] == y)],
                     False: lambda x, y: mut_dat.loc[mut_dat[x] == y]}
            method = cases[all_except]

        return method(variant, item)


    def _single_entity_mutated(self, mut_dat, output, variant, item, translocations, fusions, all_except):
        """Retrieves mutations related to single entities.
        Single entities are Gene and CellLine Classes.
        """

        out_dict = {"names": lambda x: list(set(x[self.by[self.version]])), #functions for returning specific data types
                    "dataframe": lambda x: x,
                    "dict": lambda x: dict(zip(x[self.by[self.version]], x[variant]))}

        return out_dict[output](mut_dat)


    def _multiple_entity_mutated(self, mut_dat, output, variant, item, translocations, fusions, all_except):
        """Retrieves mutatioens related to multiple entities.
        Multiple entities are Organelle
        Cancer, CellLineCluster, and GeneCluster classes.
        """

        if self.version == "canc":
            variant = "gene"
        else:
            variant = "DepMap_ID"

        out_dict = {"names": lambda x: list(set(x[self.by[self.version]])), #functions for returning specific data types
                    "dataframe": lambda x: x}

        if output == "dict":
            out = {k:mut_dat[self.by[self.version]].loc[v].unique() for k,v in mut_dat.groupby(variant).groups.items()}
        else:
            out = out_dict[output](mut_dat)

        return out

