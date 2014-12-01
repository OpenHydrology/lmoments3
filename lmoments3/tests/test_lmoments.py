import unittest
from numpy.testing import assert_almost_equal
import lmoments3 as lm
from lmoments3 import stats
from lmoments3 import distr
from lmoments3.tests import DistributionTestCase


class TestLmoments(unittest.TestCase):
    def test_samlmu(self):
        testdata = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
        expected = [3.23636364, 1.14181818, 0.27388535, 0.02335456, -0.04246285]
        result = lm.lmom_ratios(testdata)
        assert_almost_equal(result, expected)

    def test_aic(self):
        data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.4, 3.5, 1.4, 6.5, 1.2, 6.8, 5.4, 3.4]
        gamfit = distr.gam.lmom_fit(data)
        test1 = stats.AIC(data, 'gam', gamfit)
        # TODO: assert something


class TestExp(DistributionTestCase):
    dist = 'exp'
    paras = {'loc': 0.9527273, 'scale': 2.2836364}
    correct_fit = [0.9527273, 2.2836364]
    correct_qua = [1.462306, 2.535623, 4.628098]
    correct_lmr = [3.2363636, 1.1418182, 0.3333333, 0.1666667]
    correct_cdf = [0.3678311, 0.8300571, 0.9543151]
    correct_pdf = [0.11530621, 0.07441766, 0.04802853, 0.03099720]
    correct_lmom = [3.2363636, 1.1418182, 0.3806061, 0.1903030, 0.1141818]
    correct_nlogl = 20.0834598


class TestGam(DistributionTestCase):
    dist = 'gam'
    nmom = 4
    # The location parameter is not used for the gamma distribution, so len(paras) should actually return 2, for example
    # for calculating degrees of freedom etc.
    # TODO: find a nice way of dealing with this issue
    paras = {'a': 2.295206, 'loc': 0, 'scale': 1.410054}
    correct_fit = [2.295206, 1.410054]
    correct_qua = [1.447838, 2.780422, 4.766705]
    correct_lmr = [3.2363636, 1.1418181, 0.2186287, 0.1387734]
    correct_cdf = [0.3278764, 0.8222726, 0.9653452]
    correct_pdf = [0.13789672, 0.09058866, 0.05644576, 0.03391116]
    correct_lmom = [3.2363636, 1.1418181, 0.2496342, 0.1584540]
    correct_nlogl = 21.2839951


class TestGev(DistributionTestCase):
    dist = 'gev'
    paras = {'loc': 2.1792884, 'scale': 1.3956404, 'c': -0.1555609}
    correct_fit = [2.1792884, 1.3956404, -0.1555609]
    correct_qua = [1.539112, 2.705672, 4.537048]
    correct_lmr = [3.2363636, 1.1418182, 0.2738854, 0.1998461]
    correct_cdf = [0.3202800, 0.8415637, 0.9606184]
    correct_pdf = [0.13388363, 0.07913309, 0.04637529, 0.02757391]
    correct_lmom = [3.2363636, 1.1418182, 0.3127273, 0.2281879, 0.1209980]
    correct_nlogl = 21.2807575


class TestGlo(DistributionTestCase):
    dist = 'glo'
    paras = {'loc': 2.7406580, 'scale': 1.0060517, 'k': -0.2738854}
    correct_fit = [2.7406580, 1.0060517, -0.2738854]
    correct_qua = [1.580189, 2.740658, 4.437061]
    correct_lmr = [3.2363636, 1.1418182, 0.2738854, 0.2291777]
    correct_cdf = [0.3052960, 0.8519915, 0.9624759]
    correct_pdf = [0.14033225, 0.07760825, 0.04294245, 0.02463028]
    correct_lmom = [3.2363636, 1.1418182, 0.3127273, 0.2616792, 0.1439873]
    correct_nlogl = 21.6538832


class TestGno(DistributionTestCase):
    dist = 'gno'
    paras = {'loc': 2.6888917, 'scale': 1.7664322, 'k': -0.5707506}
    correct_fit = [2.6888917, 1.7664322, -0.5707506]
    correct_qua = [1.508372, 2.688892, 4.597378]
    correct_lmr = [3.2363636, 1.1418182, 0.2738848, 0.1818274]
    correct_cdf = [0.3295539, 0.8357710, 0.9599970]
    correct_pdf = [0.13099439, 0.08020770, 0.04842820, 0.02933547]
    correct_lmom = [3.2363636, 1.1418182, 0.3127266, 0.2076140, 0.1104614]
    correct_nlogl = 21.1009750


class TestGpa(DistributionTestCase):
    dist = 'gpa'
    paras = {'loc': 0.7928727, 'scale': 2.7855796, 'c': -0.1400000}  # Note: `c` parameter needs to be negated!
    correct_fit = [0.7928727, 2.7855796, 0.1400000]
    correct_qua = [1.404848, 2.632964, 4.806899]
    correct_lmr = [3.2363636, 1.1418182, 0.2738854, 0.1230499]
    correct_cdf = [0.3604888, 0.8167330, 0.9597484]
    correct_pdf = [0.12194724, 0.08343282, 0.05567275, 0.03610412]
    correct_lmom = [3.23636364, 1.14181818, 0.31272727, 0.14050066, 0.07817741]
    correct_nlogl = 20.5492690


class TestGum(DistributionTestCase):
    dist = 'gum'
    paras = {'loc': 2.285519, 'scale': 1.647295}
    correct_fit = [2.285519, 1.647295]
    correct_qua = [1.501596, 2.889274, 4.756363]
    correct_lmr = [3.236363, 1.141818, 0.169925, 0.150375]
    correct_cdf = [0.3044484, 0.8249232, 0.9693322]
    correct_pdf = [0.15060460, 0.09638151, 0.05733088, 0.03276992]
    correct_lmom = [3.2363635, 1.1418182, 0.1940235, 0.1717009, 0.0637914]
    correct_nlogl = 21.7261232


class TestKap(DistributionTestCase):
    dist = 'kap'
    paras = {'loc': -9.0633543, 'scale': 17.0127900, 'k': 0.9719618, 'h': 2.4727933}
    correct_fit = [-9.0633543, 17.0127900, 0.9719618, 2.4727933]
    correct_qua = [1.311688, 2.454434, 5.286237]
    correct_lmr = [3.23636364, 1.14181818, 0.27388545, 0.02335466]
    correct_cdf = [0.4185230, 0.7772538, 0.9769973]
    correct_pdf = [0.09794121, 0.08128701, 0.07022161, 0.06197593]
    # This cannot be correct if lmr is correct
    #correct_lmom = [3.236364, 1.141818, 0.09662925, 0.008239735, 0.00005919404]
    correct_lmom = [3.236364, 1.141818, 0.31272738604748096, 0.0266667753757188, 1.91576195e-04]

    correct_nlogl = 18.2664498


class TestNor(DistributionTestCase):
    dist = 'nor'
    paras = {'loc': 3.236364, 'scale': 2.023820}
    correct_fit = [3.236364, 2.023820]
    correct_qua = [1.533074, 3.236364, 4.939654]
    correct_lmr = [3.2363636, 1.1418182, 0.0000000, 0.1226017]
    correct_cdf = [0.2706309, 0.8082428, 0.9907083]
    correct_pdf = [0.18357864, 0.13484509, 0.07759170, 0.03497539]
    correct_lmom = [3.2363636, 1.1418182, 0.0000000, 0.1399889, 0.0000000]
    correct_nlogl = 22.6955677


class TestPe3(DistributionTestCase):
    dist = 'pe3'
    nmom = 4
    paras = {'loc': 3.236364, 'scale': 2.199489, 'skew': 1.646184}
    correct_fit = [3.236364, 2.199489, 1.646184]
    correct_qua = [1.447672, 2.663015, 4.705896]
    correct_lmr = [3.2363636, 1.1418182, 0.2738845, 0.1498865]
    correct_cdf = [0.3462110, 0.8258929, 0.9597978]
    correct_pdf = [0.12681904, 0.08243435, 0.05226950, 0.03260397]
    correct_lmom = [3.2363636, 1.1418182, 0.3127263, 0.1711432]
    correct_nlogl = 20.8203581


class TestWak(DistributionTestCase):
    dist = 'wak'
    paras = {'loc': 0.7928727, 'scale': 2.7855796, 'beta': 0.14, 'gamma': 0, 'delta': 0}
    correct_fit = [0.7928727, 2.7855796, 0.1400000, 0.0000000, 0.0000000]
    correct_qua = [1.404848, 2.632964, 4.806899]
    correct_lmr = [3.2363636, 1.1418182, 0.2738854, 0.1230499]
    correct_cdf = [0.3604888, 0.8167330, 0.9597484]
    correct_pdf = [0.12194724, 0.08343282, 0.05567275, 0.03610412]
    correct_lmom = [3.23636364, 1.14181818, 0.31272727, 0.14050066, 0.07817741]
    correct_nlogl = 20.549269


class TestWei(DistributionTestCase):
    dist = 'wei'
    paras = {'loc': 0.6740393, 'scale': 2.7087887, 'c': 1.1750218}
    correct_fit = [0.6740393, 2.7087887, 1.1750218]
    correct_qua = [1.429808, 2.656981, 4.735337]
    correct_lmr = [3.2363636, 1.1418182, 0.2738853, 0.1413359]
    correct_cdf = [0.3507727, 0.8233116, 0.9600031]
    # This appears to be incorrect
    #correct_pdf = [0.07149587, 0.04550752, 0.02836919, 0.01738135]
    # New values taken from R, using `dweibull(x-0.6740393, shape=1.1750218, scale=2.7087887)`
    correct_pdf = [0.12592833, 0.08318826, 0.05339457, 0.03351231]
    correct_lmom = [3.2363636, 1.1418182, 0.31272723, 0.16137985, -0.09159867]
    correct_nlogl = 20.749053
