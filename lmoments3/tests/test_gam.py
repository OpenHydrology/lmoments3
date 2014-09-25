import lmoments3 as lm
from lmoments3.tests import AbstractDistributionTestCase


class TestGam(AbstractDistributionTestCase):
    correct_fit = [2.295206, 1.410054]
    correct_qua = [1.447838, 2.780422, 4.766705]
    correct_lmr = [3.2363636, 1.1418181, 0.2186287, 0.1387734]
    correct_cdf = [0.3278764, 0.8222726, 0.9653452]
    correct_pdf = [0.13789672, 0.09058866, 0.05644576, 0.03391116]
    correct_lmom = [3.2363636, 1.1418181, 0.2496342, 0.1584540]

    def test_fit(self):
        fit = lm.pelgam(self.lmu)
        self.assertAlmostEqual(fit, self.correct_fit)

    def test_qua(self):
        qua = [lm.quagam(0.2, self.correct_fit),
               lm.quagam(0.5, self.correct_fit),
               lm.quagam(0.8, self.correct_fit)]
        qua2 = lm.quagam([0.2, 0.5, 0.8], self.correct_fit)
        self.assertAlmostEqual(qua, self.correct_qua)
        self.assertAlmostEqual(qua2, self.correct_qua)

    def test_lmr(self):
        lmr = lm.lmrgam(self.correct_fit, 4)
        self.assertAlmostEqual(lmr, self.correct_lmr)

    def test_cdf(self):
        cdf = [lm.cdfgam(2,self.correct_fit),
                  lm.cdfgam(5,self.correct_fit),
                  lm.cdfgam(8,self.correct_fit)]
        cdf2 = lm.cdfgam([2,5,8],self.correct_fit)
        self.assertAlmostEqual(cdf, self.correct_cdf)
        self.assertAlmostEqual(cdf2, self.correct_cdf)

    def test_pdf(self):
        pdf = [lm.pdfgam(4,self.correct_fit),
               lm.pdfgam(5,self.correct_fit),
               lm.pdfgam(6,self.correct_fit),
               lm.pdfgam(7,self.correct_fit)]
        pdf2 = lm.pdfgam([4,5,6,7],self.correct_fit)
        self.assertAlmostEqual(self.correct_pdf, pdf)
        self.assertAlmostEqual(self.correct_pdf, pdf2)

    def test_lmom(self):
        lmom = lm.lmomgam(self.correct_fit)
        self.assertAlmostEqual(lmom, self.correct_lmom)

    def test_rand(self):
        rand = lm.randgam(10,self.correct_fit)
        self.assertEqual(len(rand), 10)
