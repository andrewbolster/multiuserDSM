'''
Created on 12 Mar 2011

@author: bolster
'''
import numpy
import scipy.misc
import multiprocessing

from pycuda.tools import make_default_context
import pycuda.autoinit
import pycuda.driver as drv
import pycuda.gpuarray as garray
import pycuda.tools as pytools
from pycuda.compiler import SourceModule

kernels=SourceModule("""
//Calc PSD for all channels in parallel (probably fucks memory accesses)
__global__ void precompute_calc_psd(int[] *bitload, float gamma){
   
""")