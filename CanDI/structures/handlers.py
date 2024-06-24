import operator
import pandas as pd
import numpy as np
from collections.abc import Iterable
import six


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


        if isinstance(item, Iterable) and not isinstance(item, six.string_types):
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

