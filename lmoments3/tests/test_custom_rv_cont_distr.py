import unittest
from lmoments3 import distr
from numpy.testing import assert_almost_equal


class TestKap(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        correct_fit = {'loc': -9.0633543,
                       'scale': 17.0127900,
                       'k': 0.9719618,
                       'h': 2.4727933}
        cls.distribution = distr.kap(**correct_fit)

    def test_qua_twenty(self):
        result = self.distribution.ppf(0.2)
        expected = 1.311688
        self.assertAlmostEqual(result, expected, places=6)

    def test_qua_list(self):
        result = self.distribution.ppf([0.2, 0.5, 0.8])
        expected = [1.311688, 2.454434, 5.286237]
        assert_almost_equal(result, expected, decimal=6)

    def test_pdf_list(self):
        result = self.distribution.pdf([4, 5, 6, 7])
        expected = [0.09794121, 0.08128701, 0.07022161, 0.06197593]
        assert_almost_equal(result, expected, decimal=6)

    def test_cdf_list(self):
        result = self.distribution.cdf([2, 5, 8])
        expected = [0.4185230, 0.7772538, 0.9769973]
        assert_almost_equal(result, expected, decimal=6)


class TestWak(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        correct_fit = {'loc': 0.7928727,
                       'scale': 2.7855796,
                       'beta': 0.14,
                       'gamma': 0,
                       'delta': 0}
        cls.distribution = distr.wak(**correct_fit)

    def test_qua_twenty(self):
        result = self.distribution.ppf(0.2)
        expected = 1.404848
        self.assertAlmostEqual(result, expected, places=6)

    def test_qua_list(self):
        result = self.distribution.ppf([0.2, 0.5, 0.8])
        expected = [1.404848, 2.632964, 4.806899]
        assert_almost_equal(result, expected, decimal=6)

    def test_pdf_list(self):
        result = self.distribution.pdf([4, 5, 6, 7])
        expected = [0.12194724, 0.08343282, 0.05567275, 0.03610412]
        assert_almost_equal(result, expected, decimal=6)

    def test_cdf_list(self):
        result = self.distribution.cdf([2, 5, 8])
        expected = [0.3604888, 0.8167330, 0.9597484]
        assert_almost_equal(result, expected, decimal=6)
