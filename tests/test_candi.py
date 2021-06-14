import unittest
import pandas as pd
import numpy as np
from CanDI.candi import Entity 

class testEntity(unittest.TestCase):

    def setUp(self):

        print("setUp")
        self.gene = Entity("gene")
        self.org = Entity("org")
        self.line = Entity("line")
        self.canc = Entity("canc")
        self.matrix = pd.DataFrame(np.random.rand(17000, 800))
        self.threshold = 1
        self.return_lines = False
        self.xidx = str(np.random.randint(0, 17000))
        self.yidx = str(np.random.randint(0, 800))

    def tearDown(self):
        print("tearDown")

    def test_gene_filters(self):
        print("test_dependency_filter")

        values = self.matrix.loc[self.xidx, ]
        over = self.gene._dependency_filter(values,  'bool', "over", self.threshold, self.return_lines)
        under = self.gene._dependency_filter(values,  'bool', "under", self.threshold, self.return_lines)

        self.assertIsInstance(over, list)
        self.assertIsInstance(under, list)

        over = self.gene._dependency_filter(values, "values", "over", self.threshold, self.return_lines)
        under = self.gene._dependency_filter(values, "values", "over", self.threshold, self.return_lines)

        print(over)

