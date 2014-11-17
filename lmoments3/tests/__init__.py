import unittest
import lmoments3 as lm
from lmoments3 import distr
from numpy.testing import assert_almost_equal


class DistributionTestCase(unittest.TestCase):
    dist = None
    distr_f = None
    inputs_qua = [0.2, 0.5, 0.8]
    inputs_cdf = [2, 5, 8]
    inputs_pdf = [4, 5, 6, 7]

    @classmethod
    def setUpClass(cls):
        testdata = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
        cls.lmu = lm.samlmu(testdata)
        try:
            cls.distr_f = getattr(distr, cls.dist)
        except TypeError:
            cls.distr_f = None
        super(DistributionTestCase, cls).setUpClass()

    def assertAlmostEqual(self, first, second, places=6):
        return assert_almost_equal(first, second, decimal=places)

    def test_fit(self):
        if self.dist:
            f = getattr(lm, 'pel' + self.dist)
            fit = f(self.lmu)
            self.assertAlmostEqual(fit, self.correct_fit)

    def test_qua(self):
        if self.dist:
            f = getattr(lm, 'qua' + self.dist)
            qua = [f(self.inputs_qua[0], self.correct_fit),
                   f(self.inputs_qua[1], self.correct_fit),
                   f(self.inputs_qua[2], self.correct_fit)]
            qua2 = f(self.inputs_qua, self.correct_fit)
            self.assertAlmostEqual(qua, self.correct_qua)
            self.assertAlmostEqual(qua2, self.correct_qua)

    def test_lmr(self):
        if self.dist:
            f = getattr(lm, 'lmr' + self.dist)
            lmr = f(self.correct_fit, 4)
            self.assertAlmostEqual(lmr, self.correct_lmr)

    def test_cdf(self):
        if self.dist:
            f = getattr(lm, 'cdf' + self.dist)
            cdf = [f(self.inputs_cdf[0], self.correct_fit),
                   f(self.inputs_cdf[1], self.correct_fit),
                   f(self.inputs_cdf[2], self.correct_fit)]
            cdf2 = f(self.inputs_cdf, self.correct_fit)
            self.assertAlmostEqual(cdf, self.correct_cdf)
            self.assertAlmostEqual(cdf2, self.correct_cdf)

    def test_pdf(self):
        if self.dist:
            f = getattr(lm, 'pdf' + self.dist)
            pdf = [f(self.inputs_pdf[0], self.correct_fit),
                   f(self.inputs_pdf[1], self.correct_fit),
                   f(self.inputs_pdf[2], self.correct_fit),
                   f(self.inputs_pdf[3], self.correct_fit)]
            pdf2 = f(self.inputs_pdf, self.correct_fit)
            self.assertAlmostEqual(self.correct_pdf, pdf)
            self.assertAlmostEqual(self.correct_pdf, pdf2)

    def test_lmom(self):
        if self.dist:
            f = getattr(lm, 'lmom' + self.dist)
            lmom = f(self.correct_fit)
            self.assertAlmostEqual(lmom, self.correct_lmom)

    def test_rand(self):
        if self.dist:
            f = getattr(lm, 'rand' + self.dist)
            rand = f(10, self.correct_fit)
            self.assertEqual(len(rand), 10)


