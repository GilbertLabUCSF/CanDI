import unittest
import pandas as pd
import numpy as np
from CanDI.structures.entity import Entity 


class testEntity(unittest.TestCase):

    def setUp(self):

        print("setUp")
        #Instantiate all types of entities
        self.gene = Entity("gene")
        self.org = Entity("org")
        self.line = Entity("line")
        self.canc = Entity("canc")
        #Mock data to test on
        self.matrix = pd.DataFrame(np.random.rand(1700, 800))
        #Random single entity indeces
        self.xidx = np.random.randint(0, 1700)
        self.yidx = np.random.randint(0, 800)
        #Random multiple entity indeces
        self.xidxs = np.random.randint(0, 1700, size = np.random.randint(0, 100))
        self.yidxs = np.random.randint(0, 800, size = np.random.randint(0, 50))

    def tearDown(self):
        print("tearDown")

    def test_gene_filters(self):
        print("test_gene_class_filters")

        values = self.matrix.loc[self.xidx, ]
        over = self.gene._dependency_filter(values,  'bool', "over", 1, False)
        under = self.gene._dependency_filter(values,  'bool', "under", 1)

        self.assertIsInstance(over, list)
        self.assertIsInstance(under, list)

        over = self.gene._dependency_filter(values, "values", "over", 1, False)
        under = self.gene._dependency_filter(values, "values", "over", 1, False)

        self.assertIsInstance(over, pd.core.series.Series)
        self.assertIsInstance(under, pd.core.series.Series)

    def test_line_filters(self):
        print("test_line_class_filters")

        values = self.matrix.loc[: , self.yidx]
        over = self.line._dependency_filter(values, "bool", "over", 1, False)
        under = self.line._dependency_filter(values, "bool", "under", 1, False)

        self.assertIsInstance(over, list)
        self.assertIsInstance(under, list)

        over = self.line._dependency_filter(values, "values", "over", 1, False)
        under = self.line._dependency_filter(values, "values", "under", 1, False)

        self.assertIsInstance(over, pd.core.series.Series)
        self.assertIsInstance(under, pd.core.series.Series)


    def test_org_filters(self):
        print("test_organelle_class_filters")

        values = self.matrix.loc[self.xidxs, ]
        over = self.org._dependency_filter(values, "bool", "over", 0.5, False)
        under = self.org._dependency_filter(values, "bool", "under", 0.5, False)

        self.assertIsInstance(over, list)
        self.assertIsInstance(under, list)

        over = self.org._dependency_filter(values, "values", "over", 0.5, False)
        under = self.org._dependency_filter(values, "values", "under", 0.5,  False)

        self.assertIsInstance(over, pd.core.frame.DataFrame)
        self.assertIsInstance(under, pd.core.frame.DataFrame)


    def test_canc_filters(self):
        print("test_organelle_class_filters")

        values = self.matrix.loc[:, self.yidxs]
        over = self.canc._dependency_filter(values, "bool", "over", 0.5, False)
        under = self.canc._dependency_filter(values, "bool", "under", 0.5, False)

        self.assertIsInstance(over, list)
        self.assertIsInstance(under, list)

        over = self.canc._dependency_filter(values, "values", "over", 0.5, False)
        under = self.canc._dependency_filter(values, "values", "under", 0.5,  False)

        self.assertIsInstance(over, pd.core.frame.DataFrame)
        self.assertIsInstance(under, pd.core.frame.DataFrame)



