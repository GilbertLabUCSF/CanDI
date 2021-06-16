import pandas as pd
from CanDI.structures import handlers

class Entity(object):

    """
    Entity is the parent class for all CanDI classes.
    It's purpose is to provide a universal connection to the data for it's child classes.
    It's methods are idiomatic functions for the various handlers that apply manipulations to the datasets.

    The following functions handle most common biologically relevant queries of candi objects.
    They automatically call the filtering objects that are defined during instantiation.
    - :func:`dependent <candi.Entity.dependent>`
    - :func:`non_dependent <candi.Entity.non_depedent>`
    - :func:`dependency_of <candi.Entity.dependency_of>`
    - :func:`expressed <candi.Entity.expressed>`
    - :func:`unexpressed <candi.Entity.unexpressed>`
    - :func:`expression_of <candi.Entity.expression_of>`
    - :func:`essential <candi.Entity.essential>`
    - :func:`non_essential <candi.Entity.non_essential>`
    - :func:`effect_of <candi.Entity.effect_of>`
    - :func:`duplication <candi.Entity.duplication>`
    - :func:`deletion <candi.Entity.deletion>`
    - :func:`cn_normal <candi.Entity.cn_normal>`
    - :func:`mutated <candi.Entity.mutated>`
    """

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
            item: str or list, optional
                name or list of names of items to query expression. Can be gene symbols or DepMap_IDs
            style: str, optional
                style argument defines the type of output desired. Options are "bool" or "values"
                "bool" returns the gene/cellline name(s). "values" returns a slice from the expression dataset
            threshold: float between 0 and 1, optional
                This is the percentage of gene(s)/cellline(s) that need to pass expressed filter.
                Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
            return_lines: bool, optional
                not implemented
        Returns:
            bool
                Returns bool if item input is str
            list
                Returns list if item arg is list
            numpy.float64
                Returns if style arg is 'values' and filter result is 1-d with length of 1.
            pandas.core.series.Series
                Returns if style arg is 'values' and filter result is 1-d with length > 1.
            pands.core.frame.DataFrame
                Returns if style arg is 'values' and filter result is 2-d.
        """
        values = self._subset_handler(item, self.expression)
        return self._expression_filter(values, style, "over", threshold, return_lines)

    def unexpressed(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Unexpressed function returns genes/cellline(s) that are below a certain expression filter.

        Args:
            item: str or list, optional
                name or list of names of items to query expression. Can be gene symbols or DepMap_IDs
            style: str, optional
                style argument defines the type of output desired. Options are "bool" or "values"
                "bool" returns the gene/cellline name(s). "values" returns a slice from the expression dataset
            threshold: float between 0 and 1, optional
                This is the percentage of gene(s)/cellline(s) that need to pass unexpressed filter.
                Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
            return_lines: bool, optional
                not implemented
        Returns:
            bool
                Returns bool if item input is str
            list
                Returns list if item arg is list
            numpy.float64
                Returns if style arg is 'values' and filter result is 1-d with length of 1.
            pandas.core.series.Series
                Returns if style arg is 'values' and filter result is 1-d with length > 1.
            pands.core.frame.DataFrame
                Returns if style arg is 'values' and filter result is 2-d.
        """
        values = self._subset_handler(item, self.expression)

        return self._expression_filter(values, style, "under", threshold, return_lines)

    def expression_of(self, items):
        """It returns the expression value in (TPM) of a specific gene(s)/cellline(s).

        Args:
            items: str or list
                name of gene(s)/cellline(s) to return expression values for
        Returns:
            numpy.float64
                Returns if items arg is type str
            pandas.core.series.Series
                Returns if items arg is type list
        """
        return self._subset_handler(items, self.expression)
        # return self.expression.reindex(genes, axis=axis)
    
    def effect_of(self, items):
        """Returns the gene effect of specific gene(s)/cellline(s).
        
        Args:
            items: str or list
                name of gene(s)/cellline(s) to gene effect values for
        Returns:
            numpy.float64
                Returns if items arg is type str
            pandas.core.series.Series
                Returns if items arg is type list
        """
        return self._subset_handler(items, self.gene_effect)

    def essential(self, item=None, style="bool", threshold=1.0, return_lines=False):
        """Returns genes/cellline(s) who's gene effect is less than -1.

        Args:
            item: str or list, optional
                name or list of names of items to query gene effect. Can be gene symbols or DepMap_IDs
            style: str, optional
                style argument defines the type of output desired. Options are "bool" or "values"
                "bool" returns the gene/cellline name(s). "values" returns a slice from the gene effect dataset.
            threshold: float between 0 and 1, optional
                This is the percentage of gene(s)/cellline(s) that need to pass essential filter.
                Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
            return_lines: bool, optional
                not implemented
        Returns:
            bool
                Returns bool if item input is str
            list
                Returns list if item arg is list
            numpy.float64
                Returns if style arg is 'values' and filter result is 1-d with length of 1.
            pandas.core.series.Series
                Returns if style arg is 'values' and filter result is 1-d with length > 1.
            pands.core.frame.DataFrame
                Returns if style arg is 'values' and filter result is 2-d.
        """
        values = self._subset_handler(item, self.gene_effect)
        return self._essentiality_filter(values, style, "under", threshold, return_lines)

    def non_essential(self, item=None, style="bool", threshold=1.0, return_lines=False):
        """Returns genes/cellline(s) who's gene effect is greater than -1.
        
        Args:
            item: str or list, optional
                name or list of names of items to query gene effect. Can be gene symbols or DepMap_IDs
            style: str, optional
                style argument defines the type of output desired. Options are "bool" or "values"
                "bool" returns the gene/cellline name(s). "values" returns a slice from the gene effect dataset.
            threshold: float between 0 and 1, optional
                This is the percentage of gene(s)/cellline(s) that need to pass non_essential filter.
                Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
            return_lines: bool, optional
                not implemented
        Returns:
            bool
                Returns bool if item input is str
            list
                Returns list if item arg is list
            pandas.core.series.Series
                Returns if style arg is 'values'.
            pands.core.frame.DataFrame
                Returns if style arg is 'values'. Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
        """
        values = self._subset_handler(item, self.gene_effect)
        return self._essentiality_filter(values, style, "over", threshold, return_lines)

    def dependency_of(self, items):
        """Returns gene dependency of given items
        
        Args:
            items: str or list
                name of gene(s)/cellline(s) to gene dependency values for
        Returns:
            numpy.float64
                Returns if items arg is type str
            pandas.core.series.Series
                Returns if items arg is type list
        """
        return self._subset_handler(items, self.gene_dependency)

    def dependent(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns genes/cellline(s) whose gene dependency is greater than 0.5
        
        Args:
            item: str or list, optional
                name or list of names of items to query gene dependency. Can be gene symbols or DepMap_IDs
            style: str, optional
                style argument defines the type of output desired. Options are "bool" or "values"
                "bool" returns the gene/cellline name(s). "values" returns a slice from the gene dependency dataset.
            threshold: float between 0 and 1, optional
                This is the percentage of gene(s)/cellline(s) that need to pass dependent filter.
                Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
            return_lines: bool, optional
                not implemented
        Returns:
            bool
                Returns bool if item arg is type str
            list
                Returns list if item arg is type list
            numpy.float64
                Returns if style arg is 'values' and filter result is 1-d with length of 1.
            pandas.core.series.Series
                Returns if style arg is 'values' and filter result is 1-d with length > 1.
            pands.core.frame.DataFrame
                Returns if style arg is 'values' and filter result is 2-d.
        """
        values = self._subset_handler(item, self.gene_dependency)
        return self._dependency_filter(values, style, "over", threshold, return_lines)

    def non_dependent(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns genes/celline(s) whose gene dependency is less than 0.5
        
        Args:
            item: str or list, optional
                name or list of names of items to query gene dependency. Can be gene symbols or DepMap_IDs
            style: str, optional
                style argument defines the type of output desired. Options are "bool" or "values"
                "bool" returns the gene/cellline name(s). "values" returns a slice from the gene dependency dataset.
            threshold: float between 0 and 1, optional
                This is the percentage of gene(s)/cellline(s) that need to pass non_dependent filter.
                Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
            return_lines: bool, optional
                not implemented
        Returns:
            bool
                Returns bool if item arg is type str
            list
                Returns list if item arg is type list
            numpy.float64
                Returns if style arg is 'values' and filter result is 1-d with length of 1.
            pandas.core.series.Series
                Returns if style arg is 'values' and filter result is 1-d with length > 1.
            pands.core.frame.DataFrame
                Returns if style arg is 'values' and filter result is 2-d.
        """
        values = self._subset_handler(item, self.gene_dependency)
        return self._dependency_filter(values, style, "under", threshold, return_lines)

    def duplication(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns gene(s)/cellline(s) with copy number above specific threshold.

        Args:
            item: str or list, optional
                name or list of names of items to query gene copy number. Can be gene symbols or DepMap_IDs
            style: str, optional
                style argument defines the type of output desired. Options are "bool" or "values"
                "bool" returns the gene/cellline name(s). "values" returns a slice from the gene copy number dataset.
            threshold: float between 0 and 1, optional
                This is the percentage of gene(s)/cellline(s) that need to pass duplication filter.
                Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
            return_lines: bool, optional
                not implemented
        Returns:
            bool
                Returns bool if item arg is type str
            list
                Returns list if item arg is type list
            numpy.float64
                Returns if style arg is 'values' and filter result is 1-d with length of 1.
            pandas.core.series.Series
                Returns if style arg is 'values' and filter result is 1-d with length > 1.
            pands.core.frame.DataFrame
                Returns if style arg is 'values' and filter result is 2-d.
        """
        values = self._subset_handler(item, self.gene_cn)
        return self._copy_number_dup(values, style, "over", threshold, return_lines)

    def deletion(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns gene(s)/cellline(s) with copy number below specific threshold.

        Args:
            item: str or list, optional
                name or list of names of items to query gene copy number. Can be gene symbols or DepMap_IDs
            style: str, optional
                style argument defines the type of output desired. Options are "bool" or "values"
                "bool" returns the gene/cellline name(s). "values" returns a slice from the gene copy number dataset.
            threshold: float between 0 and 1, optional
                This is the percentage of gene(s)/cellline(s) that need to pass deletion filter.
                Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
            return_lines: bool, optional
                not implemented
        Returns:
            bool
                Returns bool if item arg is type str
            list
                Returns list if item arg is type list
            numpy.float64
                Returns if style arg is 'values' and filter result is 1-d with length of 1.
            pandas.core.series.Series
                Returns if style arg is 'values' and filter result is 1-d with length > 1.
            pands.core.frame.DataFrame
                Returns if style arg is 'values' and filter result is 2-d.
        """
        values = self._subset_handler(item, self.gene_cn)
        return self._copy_number_del(values, style, "under", threshold, return_lines)

    def cn_normal(self, item=None, style='bool', threshold=1.0, return_lines=False):
        """Returns gene(s)/cellline(s) with normal copy number.

        Args:
            item: str or list, optional
                name or list of names of items to query gene copy number. Can be gene symbols or DepMap_IDs
            style: str, optional
                style argument defines the type of output desired. Options are "bool" or "values"
                "bool" returns the gene/cellline name(s). "values" returns a slice from the gene copy number dataset.
            threshold: float between 0 and 1, optional
                This is the percentage of gene(s)/cellline(s) that need to pass copy number normal filters.
                Only relevant for Cancer, CellLineCluster, Organelle, and GeneCluster objects.
            return_lines: bool, optional
                not implemented
        Returns:
            bool
                Returns bool if item arg is type str
            list
                Returns list if item arg is type list
            numpy.float64
                Returns if style arg is 'values' and filter result is 1-d with length of 1.
            pandas.core.series.Series
                Returns if style arg is 'values' and filter result is 1-d with length > 1.
            pands.core.frame.DataFrame
                Returns if style arg is 'values' and filter result is 2-d.
        """
        values = self._subset_handler(item, self.gene_cn)
        under = self._copy_number_dup(values, "values", "under", threshold, return_lines)
        over = self._copy_number_del(under, style, "over", threshold, return_lines)
        return over

    def mutated(self, subset=None, output="names", variant=None, item=None, translocations=False, fusions=False,
                all_except=False):
        """Returns gene(s)/cellline(s) that are mutated. User can specify the type of mutation.

        Args:
            subset: str or list, optional
                name or list of names to query for mutations
            output: str
                desired datatype for output. Can be 'names', 'dataframe', or 'dict'
            variant: str
                Column in mutations data set for specific filtering
            item:
                Value in variant column for filtering
            translocations: bool
                not implemented
            fusions: bool
                not implemented
            all_except: bool
                output all values except those passing the given mutation filters
        Returns:
            list
                list of gene(s)/cellline(s) that pass mutation filters
            pandas.core.frame.DataFrame
                DataFrame of all mutational data for items that pass mutation filters
            dict
                Dictionary of with gene names as keys and list of depmap_ids as values.
        """
        if subset:
            mut_dat = self._get_mut_subset(self.mutations, subset)
            if mut_dat.empty: return
        else:
            mut_dat = self.mutations

        return self._mutation_handler(mut_dat, output, variant, item, translocations, fusions, all_except)
