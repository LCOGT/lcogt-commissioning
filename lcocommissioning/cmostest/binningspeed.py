import numpy as np
import torch
import sys
from astropy.io import fits
from astropy.nddata import block_reduce
import time

inputfile = sys.argv[1]

print (f"Inputfile is {inputfile}")

fittdata = fits.open(inputfile)['SCI'].data
print (f'Input Image dimensions: {fittdata.shape}')

DEVICE = torch.device("mps" if torch.has_mps else "cpu")
DEVICE = torch.device("cpu" if torch.has_cuda else DEVICE)
DEVICE='cpu'
torch.get_num_threads()
print (f'pytorch device: {DEVICE}, #Threads {torch.get_num_threads()}')

binning_factor = 4
iterations = 1000
def bin_astropy (data, factor):
    return block_reduce(data, factor, func=np.mean)

def bin_pyhtorch (data, factor):
    t = torch.tensor(np.array ([data.astype (np.float32), ]), device=DEVICE)
    m = torch.nn.AvgPool2d(factor)
    return m(t).cpu().numpy()[0]


start = time.time()
for rep in range (1):
    result = bin_astropy(fittdata,binning_factor)


print (f'astropy execution time: { (time.time() - start)/(rep+1)}')
print (f'Image dimensions: {result.shape}')
astropybin = result

start = time.time()
for rep in range (iterations):

    result = bin_pyhtorch(fittdata, binning_factor)
    pass
print (f'pyTorch execution time: { (time.time() - start)/(rep+1)}')
print (f'Image dimensions: {result.shape}')
tensorbin = result
delta = astropybin - tensorbin
fits.writeto('tensordelta.fits', delta, overwrite=True)