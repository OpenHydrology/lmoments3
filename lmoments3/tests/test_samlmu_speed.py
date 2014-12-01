import unittest
from datetime import datetime
import lmoments3 as lm
from lmoments3 import distr
import numpy as np


class TestSamlmuSpeed(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.n = 50
        kappa_distr = distr.kap(loc=-0.7, scale=2.6, k=0.9, h=1.6)
        cls.record = kappa_distr.ppf(np.random.random(cls.n * 10000))

    def setUp(self):
        self.start_time = datetime.now()

    def tearDown(self):
        duration = datetime.now() - self.start_time
        print("Test {} ran for {} s.".format(self.id(), duration.total_seconds()))

    def test_n50_nmom3(self):
        start_i = 0
        for i in range(10000):
            l1, l2, t3 = lm.lmom_ratios(self.record[start_i:start_i + self.n], nmom=3)
            t2 = l2 / l1
            start_i += self.n

    def test_n50_nmom4(self):
        start_i = 0
        for i in range(10000):
            l1, l2, t3, t4 = lm.lmom_ratios(self.record[start_i:start_i + self.n], nmom=4)
            t2 = l2 / l1
            start_i += self.n


class TestKappaSpeed(unittest.TestCase):
    def setUp(self):
        self.start_time = datetime.now()

    def tearDown(self):
        duration = datetime.now() - self.start_time
        print("Test {} ran for {} s.".format(self.id(), duration.total_seconds()))

    def test_n50(self):
        n = 50
        kappa_distr = distr.kap(loc=-0.7, scale=2.6, k=0.9, h=1.6)
        for i in range(10000):
            record = kappa_distr.ppf(np.random.random(n))