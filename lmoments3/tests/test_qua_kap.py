import unittest
import lmoments3 as lm
from numpy.testing import assert_almost_equal


class TestQuaKap(unittest.TestCase):
    correct_fit = [-9.0633543, 17.0127900, 0.9719618, 2.4727933]

    def test_qua_zero(self):
        result = lm.quakap(0, self.correct_fit)
        expected = 1.1797661891923266
        self.assertAlmostEqual(result, expected)

    def test_qua_zero_as_array(self):
        result = lm.quakap([0], self.correct_fit)
        print(result)
        expected = [1.1797661891923266]
        assert_almost_equal(result, expected, decimal=7)

    def test_qua_20(self):
        result = lm.quakap(0.2, self.correct_fit)
        expected = 1.311688
        self.assertAlmostEqual(result, expected, places=6)

    def test_qua_20_as_array(self):
        result = lm.quakap([0.2], self.correct_fit)
        print(result)
        expected = [1.311688]
        assert_almost_equal(result, expected, decimal=6)
