import unittest
import lmoments3 as lm
from numpy.testing import assert_almost_equal


class TestExp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        testdata = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
        cls.LMU = [round(moment, 6) for moment in lm.samlmu(testdata)]
        cls.correctexpfit = [0.9527273, 2.2836364]

    def test_expfit(self):
        expfit = lm.pelexp(self.LMU)
        assert_almost_equal(self.correctexpfit, expfit, decimal=6)

    def test_quaexp(self):
        expqua = [lm.quaexp(0.2, self.correctexpfit),
                  lm.quaexp(0.5, self.correctexpfit),
                  lm.quaexp(0.8, self.correctexpfit)]
        expqua2 = lm.quaexp([0.2, 0.5, 0.8], self.correctexpfit)
        correctexpqua = [1.462306, 2.535623, 4.628098]
        assert_almost_equal(correctexpqua, expqua, decimal=6)
        assert_almost_equal(correctexpqua, expqua2, decimal=6)

    def test_lmrexp(self):
        explmr = lm.lmrexp(self.correctexpfit, 4)
        correctexplmr = [3.2363636, 1.1418182, 0.3333333, 0.1666667]
        assert_almost_equal(correctexplmr, explmr, decimal=6)

    def test_expcdf(self):
        expcdf = [lm.cdfexp(2, self.correctexpfit),
                  lm.cdfexp(5, self.correctexpfit),
                  lm.cdfexp(8, self.correctexpfit)]
        expcdf2 = lm.cdfexp([2, 5, 8], self.correctexpfit)
        correctexpcdf = [0.3678311, 0.8300571, 0.9543151]
        assert_almost_equal(correctexpcdf, expcdf, decimal=6)
        assert_almost_equal(correctexpcdf, expcdf2, decimal=6)

    def test_pdfexp(self):
        exppdf = [lm.pdfexp(4, self.correctexpfit),
                  lm.pdfexp(5, self.correctexpfit),
                  lm.pdfexp(6, self.correctexpfit),
                  lm.pdfexp(7, self.correctexpfit)]
        exppdf2 = lm.pdfexp([4, 5, 6, 7], self.correctexpfit)
        correctexppdf = [0.11530621, 0.07441766, 0.04802853, 0.03099720]
        assert_almost_equal(correctexppdf, exppdf, decimal=6)
        assert_almost_equal(correctexppdf, exppdf2, decimal=6)

    def test_lmomexp(self):
        explmom = lm.lmomexp(self.correctexpfit)
        correctexplmom = [3.2363636, 1.1418182, 0.3806061, 0.1903030, 0.1141818]
        assert_almost_equal(correctexplmom, explmom)

    def test_randexp(self):
        exprand = lm.randexp(10, self.correctexpfit)
        self.assertEqual(len(exprand), 10)
