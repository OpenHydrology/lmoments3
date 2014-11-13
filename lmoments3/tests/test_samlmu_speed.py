import unittest
from datetime import datetime
import lmoments3 as lm
from lmoments3 import distr
import numpy as np


class TestSamlmuSpeed(unittest.TestCase):
    def setUp(self):
        self.start_time = datetime.now()

    def tearDown(self):
        duration = datetime.now() - self.start_time
        print("Test {} ran for {} s.".format(self.id(), duration.total_seconds()))

    def test_n50_nmom3(self):
        n = 50
        kappa_distr = distr.Kappa(loc=-0.7, scale=2.6, k=0.9, h=1.6)
        record = kappa_distr.ppf(np.random.random(n*10000))
        start_i = 0
        for i in range(10000):
            l1, l2, t3 = lm.samlmu(record[start_i:start_i + n], nmom=3)
            t2 = l2 / l1
            start_i += n


class TestKappaSpeed(unittest.TestCase):
    def setUp(self):
        self.start_time = datetime.now()

    def tearDown(self):
        duration = datetime.now() - self.start_time
        print("Test {} ran for {} s.".format(self.id(), duration.total_seconds()))

    def test_n50(self):
        n = 50
        kappa_distr = distr.Kappa(loc=-0.7, scale=2.6, k=0.9, h=1.6)
        for i in range(10000):
            record = kappa_distr.ppf(np.random.random(n))