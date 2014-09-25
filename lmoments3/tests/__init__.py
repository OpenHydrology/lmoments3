import unittest
import lmoments3 as lm
from numpy.testing import assert_almost_equal


class AbstractDistributionTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        testdata = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
        cls.lmu = lm.samlmu(testdata)

        super(AbstractDistributionTestCase, cls).setUpClass()

    def assertAlmostEqual(self, first, second, places=6):
        return assert_almost_equal(first, second, decimal=places)