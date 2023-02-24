import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.stats import rv_continuous
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import datetime
from sklearn.mixture import GaussianMixture

class Telegraph(rv_continuous):
    "negative exponential"
    def __init__(self):
        self.x0 = 420
        self.sigma = 5
        self.delta = 20

    def _pdf(self, x, x0, sigma):

        return norm (loc=self.x0, scale=self.sigma)

    def _argcheck(self, x0, sigma):
        return Truescipy.optimize


def myexp (x, A, x0,sigma):
    return A * np.exp ( -1 * (x-x0)**2 / (2*sigma**2))

def fit_function (x, A, B, x0, sigma, delta):
    return myexp (x, A, x0, sigma) + myexp (x,B,x0-delta, sigma) + myexp (x, B, x0+delta, sigma)

def summed_ln_likelihood (x, A,B,x0,sigma, delta):
    return np.sum (np.log ( fit_function(x, A, B, x0,sigma, delta)))

def findditributionparamters (xvec):
    guess = [0.5,0.25, np.mean(xvec), 5, 5 ]
    def lnprob (parameters, xvec ):
        return -1 * summed_ln_likelihood(xvec,parameters[0],parameters[1],parameters[2],parameters[3],parameters[4],)
    res = scipy.optimize.minimize(lnprob, guess, xvec,)
    print ("Likelihood fit:", res.x)

def findx0 (xvec, A, B, sigma, delta):

    def lnprob (paramters, xvec):
        #print (f"paramters: {paramters} data: {xvec}")
        x0 = paramters[0]
        return -1 * summed_ln_likelihood(xvec,A,B,x0, sigma, delta)

    start = datetime.datetime.utcnow()

    guess = [np.mean (xvec), ]
    res = scipy.optimize.minimize(lnprob, guess, xvec,)
    end = datetime.datetime.utcnow()
    return (res)

def readdistribution (file):
    return np.asarray(np.loadtxt (file))

def fitandplot (data, label):
    numbins = 20
    bins = np.linspace (np.min(data), np.max(data), numbins)
    histo1, bins1 = np.histogram(data, bins=bins, normed=True)
    binscenters = np.array([0.5 * (bins1[i] + bins1[i+1]) for i in range(len(bins1)-1)])
    popt = None
    pcov = None
    try:
        popt, pcov = curve_fit(fit_function, xdata=binscenters, ydata=histo1, p0=[0.5, 0.5, np.mean (data),5 ,40],
                               bounds = [ [0, 0, np.min(data), 0,0], [1, 1,  np.max(data), np.inf, np.inf]],
                             )
        print(popt, len(binscenters))
    except :
        print ("Fitting failed")

    plt.hist ( data, bins=bins1, density = True,   label=f"data {label}")



    # gaussian model:
    clf = GaussianMixture(n_components=3, covariance_type="full")
    clf.fit (data.reshape (-1,1))

    m1 = clf.means_[0][0]
    m2 = clf.means_[1][0]
    m3 = clf.means_[2][0]
    s1 = np.sqrt(clf.covariances_[0][0][0])
    s2 = np.sqrt(clf.covariances_[1][0][0])
    s3 = np.sqrt(clf.covariances_[2][0][0])
    w1 = clf.weights_[0]
    w2 = clf.weights_[1]
    w3 = clf.weights_[2]
    print ("GMM:", m1,m2,m3,"\n", s1,s2,s3,"\n", w1,w2,w3)

    popt = np.asarray([w1,(w2+w3)/2, m1, (s1+s2+s3)/3, np.abs(m3-m1)])
    print (popt)

    if popt is not None:
        xspace = np.linspace(np.min(data), np.max(data), 100)
        plt.plot (xspace, fit_function(xspace, *popt), label=f"Fit: {popt[0]: 4.2f} {popt[1]: 4.2f} {popt[2]: 7.1f} "
                                                         f"{popt[3]:> 4.1f} {popt[4]:> 8.1f}")

plt.figure()
data = readdistribution('exampledata/mindistr.txt')
fitandplot(data, "Best case")
findditributionparamters(data)
data = readdistribution('exampledata/maxdistr.txt')
fitandplot(data, "Worst case")
findditributionparamters(data)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
plt.savefig ("distributionfit.png", bbox_inches="tight")
plt.close()


x = range (len(data))
means = [np.mean (data[0:x]) for x in range(len(data))]
mlhres = [ findx0(data[0:x], 0.06, 0.02, 7, 37.0 ) for x in range (len(data))]
mlh = [res.x[0] for res in mlhres]
plt.figure()
plt.plot (x,data,'.')
plt.plot (x,means, label = "Mean of first N")
plt.hlines(np.mean (data), xmin=0, xmax=len(data), color='black')
plt.plot (x,mlh, label = "maximum likelihood fit")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
plt.savefig ("meanvsmle.png",bbox_inches="tight", dpi=300)

plt.figure()
plt.plot (x,data,'.')
plt.plot (x,means, label = "Mean of first N")
plt.hlines(np.mean (data), xmin=0, xmax=len(data), color='black')
plt.plot (x,mlh, label = "maximum likelihood fit")
plt.ylim([np.mean (data)-5,np.mean (data)+5])
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05))
plt.savefig ("meanvsmle_zoom.png",bbox_inches="tight", dpi=300)
plt.close()

# # timing test
# _data = data[0:25]
# start = datetime.datetime.utcnow()
# for ii in range (1000):
#     findx0(_data, 0.06, 0.02, 5, 39.3 )
# end = datetime.datetime.utcnow()
# print (f"Timing test: {(end - start) / 1000}")

plt.figure()
x = np.linspace(min(data), max (data), 1000)
y=[summed_ln_likelihood(data[0:], 0.06,0.02,_x, 5,39.3) for _x in x]
plt.plot (x,y, label="ln likelihood")
plt.savefig ("lhe.png")




