import operator
import pandas as pd
import numpy as np
from collections import MutableSequence
from CanDI import data


class Grabber:

    """"
    Grabber class handles all bulk data retrival from the onc Classes.
    """

    def __init__(self, grabber_type, key, axis):

        self.key = key
        self.gtype = getattr(self, grabber_type)
        self._isin_col = self.isin_dict[axis]
        self.axis = axis

    def __call__(self, item):

        if item not in dir(data):
            raise AttributeError("data has no attribute {}".format(item))

        dataset = getattr(data, item)

        if type(dataset) is str:

            to_load = input("{} has not been loaded. Do you want to load, y/n?> ".format(item))
            if to_load == ("y" or "Y" or "Yes"):
                dataset = data.load(item)
                print("Load Complete")
            else:
                return

        return self.gtype[item](dataset)

    @property
    def isin_dict(self):

        return {0: "gene",
                1: "DepMap_ID"}

    @property
    def gene(self):

        return {"pickles": self.get_one,
                "depmap": self.get_one,
                "transcription": self.get_one,
                "counts": self.get_one,
                "copy_number": self.get_one,
                "complexes": self.get_one,
                "locations": self.isin,
                "mutations": self.isin,
                "fusions": self.merge_two,
                "translocations": self.merge_two,
                "interactions": self.merge_two}
    @property
    def line(self):

        return {"pickles": self.get_one,
                "depmap": self.get_one,
                "transcription": self.get_one,
                "counts": self.get_one,
                "copy_number": self.get_one,
                "mutations": self.isin,
                "fusions": self.merge_two,
                "translocations": self.isin}


    @property
    def canc(self):

        return {"pickles": self.get_several,
                "depmap": self.get_several,
                "transcription": self.get_several,
                "copy_number": self.get_several,
                "mutations": self.isin,
                "fusions": self.isin,
                "translocations": self.isin,
                "counts": self.get_several}


    @property
    def org(self):

        return {"pickles": self.get_several,
                "depmap": self.get_several,
                "transcription": self.get_several,
                "copy_number": self.get_several,
                "complexes": self.get_several,
                "counts": self.get_several,
                "mutations": self.isin,
                "fusions": self.merge_two,
                "translocations": self.merge_two,
                "interactions": self.merge_two}

    def get_one(self, dataset):

        cases = {0: lambda x,y: x.loc[y],
                 1: lambda x,y: x[y]}
        try:
            series = cases[self.axis](dataset, self.key)
        except KeyError:
            return

        return series

    def get_several(self, dataset):

        func = lambda x,y: x.reindex(y, axis=self.axis)
        values = func(dataset, self.key).dropna(how="all", axis=self.axis)

        if values.empty:
            return
        else:
            return values

    def merge_two(self, dataset):

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

    def isin(self, dataset):
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






class BinaryFilter:


    def __init__(self, margin, handler):

        self.margin = margin
        self._default = self._handlers[handler]

    @property
    def _handlers(self):

        return {np.float64: self._float_handler,
                pd.Series: self._series_handler,
                pd.DataFrame: self._frame_handler}

    def __call__(self, vals, style, caller, threshold, return_lines=False):

        self._eval_args(vals, style, threshold, return_lines)

        handler = self._handlers.get(type(vals), self._default)
        return handler(vals, style, caller, threshold, return_lines)

    def _float_handler(self, values, style, caller, *args):

        if style == 'values':
            return values
        else:
            behaviors = {"over": operator.ge,
                         "under": operator.lt}

            return behaviors[caller](values, self.margin)

    def _series_handler(self, values, style, caller, *args):

        behaviors = {"over": values.ge,
                     "under": values.lt}

        evaluated = values[behaviors.get(caller)(self.margin)]

        if style == "values":
            return evaluated
        else:
            return list(evaluated.index)

    def _frame_handler(self, values, style, caller, threshold, return_lines):

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

#        try:
#            formatted = evaluated.to_dict("index").items()
#        except TypeError:
#            formatted = evaluated.unstack().to_dict("index").items()
#
#        eva_dict =  {k:list(v.keys()) for k, v in formatted}
#
#        if return_lines:
#            return eva_dict
#        else:
#            return list(eva_dict.keys())

    @staticmethod
    def _eval_args(vals=None, style="bool", threshold=1.0, return_lines=False):
        assert style in ["bool", "values"],"style must be 'bool' or 'values'"
        if return_lines:
            assert len(self._obj.shape) > 1, "return line is not available for type {}".format(self._obj)
        assert 0.0 < threshold <= 1.0, "threshold is invalid, must be between 0 and 1"
        assert vals is not None
        try:
            assert ~vals.empty
        except AttributeError:
            pass
        return











class MutationHandler(object):

    """
    Handler filtering and querying of mutation data.
    Has methods for all onc objects that allow for more specific filtering.
    Includes methods for querying translocations and fusions.
    """


    def __init__(self, version):
        self.by = {"gene":"DepMap_ID",
                   "line":"gene",
                   "canc":"DepMap_ID",
                   "org": "DepMap_ID"}
        self.version = version


    def __call__(self, mut_dat, output, variant, item, translocations, fusions, all_except):

        if variant and item:
            mut_dat = self._get_variant(mut_dat, variant, item, all_except=all_except)
        else:
            try:
                mut_dat = self._get_variant(mut_dat, "Variant_Classification", "Silent", all_except=True)
            except AssertionError:
                pass

        cases = {"gene": self._single_entity_mutated,
                "line": self._single_entity_mutated,
                "canc": self._multiple_entity_mutated,
                "org": self._multiple_entity_mutated}

        return cases[self.version](mut_dat, output, variant, item, translocations, fusions, all_except)

    @staticmethod
    def _get_variant(mut_dat, variant, item, all_except=False):

        assert item in mut_dat[variant].unique(), "{0} not found, options are: {1}".format(item, mut_dat[variant].unique())


        if isinstance(item, MutableSequence):
            method = lambda x,y: mut_dat.loc[mut_dat[x].isin(y)]

        else:
            cases = {True: lambda x, y: mut_dat.loc[~(mut_dat[x] == y)],
                     False: lambda x, y: mut_dat.loc[mut_dat[x] == y]}
            method = cases[all_except]

        return method(variant, item)


    def _single_entity_mutated(self, mut_dat, output, variant, item, translocations, fusions, all_except):

        out_dict = {"names": lambda x: list(set(x[self.by[self.version]])),
                    "dataframe": lambda x: x,
                    "dict": lambda x: dict(zip(x[self.by[self.version]], x[variant]))}

        return out_dict[output](mut_dat)



    def _multiple_entity_mutated(self, mut_dat, output, variant, item, translocations, fusions, all_except):

        if self.version == "canc":
            variant = "gene"
        else:
            variant = "DepMap_ID"

        out_dict = {"names": lambda x: list(set(x[self.by[self.version]])),
                    "dataframe": lambda x: x}

        if output == "dict":
            out = {k:mut_dat[self.by[self.version]].loc[v].unique() for k,v in mut_dat.groupby(variant).groups.items()}
        else:
            out = out_dict[output](mut_dat)

        return out




