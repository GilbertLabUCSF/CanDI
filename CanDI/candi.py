# Classes for handling data aggregations
import operator
from collections import OrderedDict, MutableSequence
import itertools as it
import pandas as pd
import numpy as np
from CanDI import data, handlers


class SubsetHandler(object):

    """
    SubsetHandler gets subsets from various datasets.
    It provides the logic for determining how to query various datasets.
    Automates finding what type of argument the user provided.
    """

    def __call__(self, arg, dat):

        if arg is None:
            return dat
        else:
            if isinstance(arg, np.ndarray):
                arg = list(arg)
            cases = {str: self.handle_string,
                     tuple: self.handle_tuple,
                     list: self.get_many}

            return cases[type(arg)](arg, dat)

    def handle_string(self, arg, dat):

        if arg in data.genes.index or arg in data.cell_lines.index:
            ids = arg

        elif arg in data.cell_lines["Primary Disease"].unique():
            ids = self.get_cancer_ids(arg)
            return self.get_many(ids, dat)

        elif arg in data.locations.location.unique():
            ids = Organelle(arg).genes
            return self.get_many(ids, dat)

        else:
            raise KeyError("{} not found in any dataset".format(arg))

        return self.get_one(ids, dat)

    def handle_tuple(self, arg, dat):

        ids = self.get_cancer_ids(arg[0], subtype=arg[1])
        return self.get_many(ids, dat)

    @staticmethod
    def get_cancer_ids(cancer, subtype=None):

        return Cancer(cancer, subtype=subtype).depmap_ids

    @staticmethod
    def get_organelle_genes(organelle, conf=6):

        return Organelle(organelle, conf=conf).genes

    @staticmethod
    def get_one(arg, df):

        try:
            vals = df.loc[arg]
        except KeyError:
            vals = df.loc[:, arg]

        return vals

    @staticmethod
    def get_many(arg, dat):

        try:
            vals = dat.reindex(arg, axis=0).dropna(how='all', axis=0)
            assert not vals.empty
            return vals

        except AssertionError:
            vals = dat.reindex(arg, axis=1).dropna(how='all', axis=1)
            assert not vals.empty
            return vals


###################################################################################################


class Entity(object):

    """
    Entity is the parent class for all CanDI classes.
    It's purpose is to provide a universal connection to the data for it's child classes.
    It's methods are idiomatic functions for the various handlers that apply manipulations to the datasets.

    The following functions handle most common biologically relevant queries of candi objects.
    They automatically call the filtering objects that are defined during instantiation.
    - :func:`dependent <candi.Entity.dependent>`
    - :func:`non_dependent <candi.Entity.non_depedent>`
    - :func:`expressed <candi.Entity.expressed>`
    - :func:`unexpressed <candi.Entity.unexpressed>`
    - :func:`expression_of <candi.Entity.expression_of>`
    - :func:`duplication <candi.Entity.duplication>`
    - :func:`deletion <candi.Entity.deletion>`
    - :func:`cn_normal <candi.Entity.cn_normal>`
    - :func:`mutated <candi.Entity.mutated>`
    """

    # - :func:`essentiality_of <candi.Entity.essentiality_of>`
    # - :func:`essential <candi.Entity.essential>`
    # - :func:`non_essential <candi.Entity.non_essential>`

    def __init__(self, obj):

        if obj == ("gene" or "org"):
            self._axis = 0
        else:
            self._axis = 1

        if obj == ("gene" or "line"):
            bi_filt = pd.Series
        else:
            bi_filt = pd.DataFrame

        """Entity functions rely on the following data handlers."""
        self._dependency_filter = handlers.BinaryFilter(0.50, bi_filt)
        self._essentiality_filter = handlers.BinaryFilter(-1.0, bi_filt)
        self._expression_filter = handlers.BinaryFilter(1.0, bi_filt)
        self._copy_number_del = handlers.BinaryFilter(0.92, bi_filt)
        self._copy_number_dup = handlers.BinaryFilter(1.07, bi_filt)
        self._mutation_handler = handlers.MutationHandler(obj)
        self._subset_handler = SubsetHandler()

    def __getattr__(self, attr):

        values = self._grabber(attr)

        setattr(self, attr, values)
        return values

    # """The following functions handle most common biologically relevant queries of candi objects.
    # They automatically call the filtering objects that are defined during instantiation.
    # """

    def expressed(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns genes/celllines that are above a certain expression filter.

        Args:
            item:
            style:
            threshold:
            return_lines:
        Returns:
            Filtered `Entity` object
        """
        values = self._subset_handler(item, self.expression)
        return self._expression_filter(values, style, "over", threshold, return_lines)

    def unexpressed(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Unexpressed function returns genes/cellline(s) that are below a certain expression filter.

        Args:
            item: str, optional
            style: str, optional
            threshold: float, optional
            return_lines: bool, optional
        Returns:
            Filtered `Entity` object
        """
        values = self._subset_handler(item, self.expression)

        return self._expression_filter(values, style, "under", threshold, return_lines)

    def expression_of(self, items):
        """It returns the expression value in (TPM) of a specific gene(s)/cellline(s).

        Args:
            items: str, list
        Returns:
            float
        """
        return self._subset_handler(items, self.expression)
        # return self.expression.reindex(genes, axis=axis)
    
    def effect_of(self, items):
        """Returns the gene effect of specific gene(s)/cellline(s).
        
        Args:
            items: str, list
        Returns:
            float
        """
        return self._subset_handler(items, self.gene_effect)

    def essential(self, item=None, style="bool", threshold=1.0, return_lines=False):
        """Returns genes/cellline(s) who's gene effect is less than -1.

        Args:
            item: str, optional
            style: str, optional
            threshold: float, optional
            return_lines: bool, optional
        Returns:
        """
        values = self._subset_handler(item, self.gene_effect)
        return self._essentiality_filter(values, style, "under", threshold, return_lines)

    def non_essential(self, item=None, style="bool", threshold=1.0, return_lines=False):
        """Returns genes/cellline(s) who's gene effect is greater than -1.
        Args:
            items:
        Returns:
        """
        values = self._subset_handler(item, self.gene_effect)
        return self._essentiality_filter(values, style, "over", threshold, return_lines)

    def dependency_of(self, items):
        """Returns gene dependency of given items
        Args:
            items:
        Returns:
        """
        return self._subset_handler(items, self.gene_dependency)

    def dependent(self, item=None, style='bool', threshold=1.0, return_lines=False):#self, item=None, style='bool'):
        """Returns genes/cellline(s) whose gene dependency is greater than 0.5
        Args:
            item: str, optional
            style: str, optional
            threshold: float, optional
            return_lines: bool, optional
        Returns:
        """
        values = self._subset_handler(item, self.gene_dependency)
        return self._dependency_filter(values, style, "over", threshold, return_lines)

    def non_dependent(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns genes/celline(s) whose gene dependency is greater than 0.5
        Args:
            item: str, optional
            style: str, optional
            threshold: float, optional
            return_lines: bool, optional
        Returns:

        """
        values = self._subset_handler(item, self.gene_dependency)
        return self._dependency_filter(values, style, "under", threshold, return_lines)

    def duplication(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns gene(s)/cellline(s) with copy number above specific threshold.

        Args:
            item:
            style:
            threshold:
            return_lines:
        Returns:
            Filtered `Entity` object
        """
        values = self._subset_handler(item, self.gene_cn)
        return self._copy_number_dup(values, style, "over", threshold, return_lines)

    def deletion(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns gene(s)/cellline(s) with copy number below specific threshold.

        Args:
            item:
            style:
            threshold:
            return_lines:
        Returns:
            Filtered `Entity` object
        """
        values = self._subset_handler(item, self.gene_cn)
        return self._copy_number_del(values, style, "under", threshold, return_lines)

    def cn_normal(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns gene(s)/cellline(s) with normal copy number.

        Args:
            item:
            style:
            threshold:
            return_lines:
        Returns:
            Filtered `Entity` object
        """
        values = self._subset_handler(item, self.gene_cn)
        under = self._copy_number_dup(values, "values", "under", threshold, return_lines)
        over = self._copy_number_del(under, style, "over", threshold, return_lines)
        return over

    def mutated(self, subset=None, output="names", variant=None, item=None, translocations=False, fusions=False,
                all_except=False):
        """Returns gene(s)/cellline(s) that are mutated. User can specify the type of mutation.

        Args:
            subset:
            output:
            variant:
            item:
            translocations:
            fusions:
            all_except:
        Returns:
            gene(s)/cellline(s) that are mutated
        """
        if subset:
            mut_dat = self._get_mut_subset(self.mutations, subset)
            if mut_dat.empty: return
        else:
            mut_dat = self.mutations

        return self._mutation_handler(mut_dat, output, variant, item, translocations, fusions, all_except)


###################################################################################################


class Gene(Entity):

    """Class used for gathering information on a single gene.
    Instantiated by gene name (preferred) or ENTREZ ID.
    Note: Not all genes have been asigned entrez ids and gene names are inconsistent across sources.
    If something doesn't show up, try alternate names.
    """

    def __init__(self, name, by="symbol"):
        super().__init__("gene")
        assert type(name) == str, "name must be string"

        if by not in ["name", "symbol", "entrez", "ensembl"]:
            raise ValueError("""by must be in ["name", "symbol", "entrez", "ensembl"] """)

        by_dict = {"name": "Approved name", "entrez": "ENTREZ ID", "ensembl": "Ensembl ID"}

        try:
            info = data.genes.loc[name]
        except (ValueError, KeyError):
            info = data.genes.loc[data.genes[by_dict[by]] == name].iloc[0]

        self.symbol = info.name
        self.name = info["Approved name"]
        self.entrez = info["ENTREZ ID"]
        self.ensembl = info["Ensembl ID"]
        self._grabber = handlers.Grabber("gene", self.symbol, self._axis)

    @property
    def get_name(self):
        return self.symbol

    def _get_mut_subset(self, mut_dat, subset):

        cases = {str: lambda x: Cancer(x).depmap_ids,
                 tuple: lambda x: Cancer(*x).depmap_ids,
                 list: lambda x: x}

        subset_ids = cases[type(subset)](subset)
        return mut_dat.loc[mut_dat.DepMap_ID.isin(subset_ids)]


###################################################################################################


class Organelle(Entity):
    """Organelle class is a group of genes defined by their subcellular location.
    It is instantiated by providing a specific organelle. During instantiation
    candi with will subset the locations dataset by the user provided organelle.
    """
    def __init__(self, organelle, min_conf=3):
        super().__init__("org")
        self.location = organelle
        self.conf = min_conf
        locs = data.locations.loc[data.locations.location == organelle]
        self.genes_and_conf = locs.loc[locs.confidence >= min_conf, :].reindex(["gene", "confidence"], axis=1)
        self.genes = list(self.genes_and_conf.gene)
        self._string_meth = lambda x, y: x[y]
        self._grabber = handlers.Grabber("org", self.genes, self._axis)

    @property
    def get_name(self):
        return self.location


###################################################################################################


class GeneCluster(Entity):
    """Functions the same as Organelle, except the genes are predetermined by the user.
    """
    def __init__(self, genes, name=None):
        super().__init__("org")
        self.genes = genes
        locs = data.locations.loc[data.locations.gene.isin(genes)]
        self._string_meth = lambda x, y: x[y]
        self.name = name


###################################################################################################


class CellLine(Entity):
    """Contains methods for gather data for a specific cell line.
    Can be instantiated by DepMap_ID (preferred) or name (in all caps).
    """
    def __init__(self, cellline):
        super().__init__("line")
        assert type(cellline) == str, "cellline must be string"

        try:
            info = data.cell_lines.loc[cellline]
        except KeyError:
            info = data.cell_lines.loc[data.cell_lines.stripped_cell_line_name == cellline].iloc[0]
        except KeyError:
            info = data.cell_lines.loc[data.cell_lines.CCLE_Name == cellline].iloc[0]
        except IndexError:
            raise ValueError("Cannot Instantiate CellLine object with {}".format(cellline))

        self.depmap_id = info.name
        self.ccle_name = info.CCLE_Name
        self.name = info["stripped_cell_line_name"]
        self.tissue = info.lineage
        self.aliases = info.Alias
        self.cosmic_id = info.COSMICID
        self.sanger_id = info["Sanger_Model_ID"]
        self.sex = info.sex
        self.source = info.source
        self.lineage = info["lineage"]
        self.subtype = info["lineage_subtype"]
        self._grabber = handlers.Grabber("line", self.depmap_id, self._axis)

    @property
    def get_name(self):
        return self.depmap_id

    def _get_mut_subset(self, mut_dat, subset):

        if type(subset) is str and subset in data.locations.location.unique():
            mut_dat_subset = mut_dat.loc[mut_dat.gene.isin(Organelle(subset).genes)]

        elif type(subset) is list:
            mut_dat_subset = mut_dat.loc[mut_dat.gene.isin(subset)]

        else:
            mut_dat_subset = mut_dat.loc[mut_dat.gene == subset]

        return mut_dat_subset


###################################################################################################


class Cancer(Entity):
    """Collects all data from cell lines of a specific disease or subtype.
    """
    def __init__(self, disease, subtype=None, gender=None, source=None, all_except=False):
        super().__init__("canc")

        if subtype is None and all_except is False:
            info = data.cell_lines.loc[data.cell_lines["primary_disease"].isin([disease])]
        elif subtype is None and all_except is True:
            info = data.cell_lines.loc[~data.cell_lines["primary_disease"].isin([disease])]
            disease = "All Except {}".format(disease)
        else:
            info = data.cell_lines.loc[data.cell_lines["lineage_subtype"] == subtype, ]
        if gender:
            info = info.loc[info.sex == gender]
        if source:
            info = info.loc[info.source == source]

        self.disease = disease
        self.depmap_ids = list(info.index)
        self.names = list(info.cell_line_name)
        self.ccle_names = list(info.CCLE_Name)
        self.subtypes = info["Subtype"].unique()
        self.sexes = info.sex.unique()
        self.sources = info.source.unique()
        self._info = info
        self._string_meth = lambda x, y: x.loc[y]
        self._grabber = handlers.Grabber("canc", self.depmap_ids, self._axis)

    @property
    def get_name(self):
        return self.disease

    def mutation_matrix(self, subset=None):
        """Returns binary nxm dataframe with genes as rows and cell lines as columns.
        Note:
            If the nth row and mth column is equal to 0 the nth gene is not mutated in the mth cell line.
            If the nth row and mth column is equal to 1 the nth gene is mutated in the mth cell line.
        Args:
            subset:
        Returns:
            binary nxm dataframe
        """
        mut_dict = self.mutated(subset, output="dict")
        in_list = []
        all_values = []
        for k, v in mut_dict.items():
            mut_dict[k] = dict(zip(v, [1] * len(v)))
            not_muts = list(set(self.depmap_ids) - set(v))
            not_muts = dict(zip(not_muts, [0] * len(not_muts)))
            mut_dict[k].update(not_muts)

        return pd.DataFrame().from_dict(mut_dict)

    def _get_mut_subset(self, mut_dat, subset, subset_type="Gene"):

        if subset_type == "Gene":
            try:
                assert isinstance(subset, MutableSequence) or isinstance(subset, np.ndarray)
            except AssertionError:
                subset = [subset]

        elif subset_type == "Organelle":
            subset = Organelle(subset).genes

        return mut_dat.loc[mut_dat.gene.isin(subset)]


###################################################################################################


class CellLineCluster(Entity):
    """Functions the same as a :class:`Cancer` object, but has predefined set of cell lines defined by user.
    """
    def __init__(self, lines, all_except=False):
        super().__init__("canc")

        assert isinstance(lines, MutableSequence), "Must be list or array-like"

        if not all_except:
            info = data.cell_lines.loc[lines]
        else:
            info = data.cell_lines[~lines]

        self.disease = info["lineage"].unique()
        self.subtypes = info["lineage_subtype"].unique()
        self.names = list(info.cell_line_name)
        self.ccle_names = list(info.CCLE_Name)
        self.depmap_ids = list(info.index)
        self.genders = info.sex.unique()
        self.sources = info.source.unique()
        self._info = info
        self._string_meth = lambda x, y: x[y]
        self._grabber = handlers.Grabber("canc", self.depmap_ids, self._axis)

    def mutation_matrix(self, subset=None):
        """Returns binary nxm dataframe with genes as rows and cell lines as columns.
        Note:
            If the nth row and mth column is equal to 0 the nth gene is not mutated in the mth cell line.
            If the nth row and mth column is equal to 1 the nth gene is mutated in the mth cell line.
        Args:
            subset:
        Returns:
            binary nxm dataframe
        """
        mut_dict = self.mutated(subset, output="dict")
        in_list = []
        all_values = []
        for k, v in mut_dict.items():
            mut_dict[k] = dict(zip(v, [1] * len(v)))
            not_muts = list(set(self.depmap_ids) - set(v))
            not_muts = dict(zip(not_muts, [0] * len(not_muts)))
            mut_dict[k].update(not_muts)

        return pd.DataFrame().from_dict(mut_dict)

    def _get_mut_subset(self, mut_dat, subset, subset_type="Gene"):

        if subset_type == "Gene":
            try:
                assert isinstance(subset, MutableSequence)
            except AssertionError:
                subset = [subset]

        elif subset_type == "Organelle":
            subset = Organelle(subset).genes

        return mut_dat.loc[mut_dat.gene.isin(subset)]

    @property
    def get_name(self):
        return self.disease
