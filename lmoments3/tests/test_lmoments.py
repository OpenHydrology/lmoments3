import unittest
import lmoments3 as lm
from numpy.testing import assert_almost_equal


class TestLmoments(unittest.TestCase):
    def test_samlmu(self):
        testdata = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
        expected = [3.23636364, 1.14181818, 0.27388535, 0.02335456, -0.04246285]
        result = lm.samlmu(testdata)
        assert_almost_equal(result, expected)