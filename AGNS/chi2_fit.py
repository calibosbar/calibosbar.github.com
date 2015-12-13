import matplotlib.pyplot as plt
import scipy
import scipy.stats
size = 10000
x = scipy.arange(size)
#y = scipy.stats.chi2.rvs(6,2,size=size)
y=np.array([1187.7377355082149,1163.7951354172355,1174.1648382707642,2101.6737019327734,2490.5802151243556,229.31686185335579,3702.9850404366912,8010.8976785589657,1426.7212688049581,731.80523340894786,633.8375833521842,668.22543662170722,461.59385929530822,1075.2450141562779,1038.3811265072629, 2069.479220105387, 1649.2137587619873,844.8898560351239,1402.7387981047777,2353.4017982222886,1152.548420343079,590.47101548905209,1849.3699518480601,1052.9411549059328])
# creating the histogram
h = plt.hist(y,normed=True, histtype='stepfilled', alpha=0.2)


dist_names = ['chi2']

for dist_name in dist_names:
    dist = getattr(scipy.stats, dist_name)
    param = dist.fit(y)
    pdf_fitted = dist.pdf(x, *param[:-2], loc=param[-2], scale=param[-1])
    plt.plot(pdf_fitted, label=dist_name)
plt.legend(loc='upper left')
plt.show()

