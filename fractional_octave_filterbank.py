# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:36:18 2015

@author: matthias
"""

import sys
from os.path import join, dirname

sys.path.append(join(dirname(__file__), "resources", "preferred_numbers"))
from preferred_numbers import preferred_number

import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as patches

from scipy.signal import spectrogram

# TODO: tested only with overlaps of 0 and 0.5
# TODO: fraction parameter not tested - i.e. not working

def get_bin_indices_for_octaves(idx_octave, L_DFT, fs, fraction=1, overlap_factor=0):
    # compute the corner frequency indices according to
    # ISO 3741:2010 (hopefully, as i don't have access to that document)
    #    f_base = 16 # the base frequency from which everything else is computed

    f_base = 16
    
    f_center = preferred_number(f_base * 2**(idx_octave/(1+2*overlap_factor)), 80)[0]
#    f_center = f_base * 2 ** (idx_octave/(1+2*overlap_factor))
 
    
    # compute lower and upper frequencies
    f_lower = f_center * 2**(-fraction/2)
    f_upper = f_center * 2**(fraction/2)
    
    f_upper = min(f_upper, fs/2)
    
    # convert to bin indices
    delta_f = fs / L_DFT
    
    idx_lower = f_lower // delta_f
    idx_center = f_center // delta_f
    idx_upper = f_upper // delta_f
    
#    print(f_lower, f_center, f_upper)
    
    return(idx_lower, idx_center, idx_upper)
    
def get_N_bands(L_DFT, fs, fraction=1, overlap_factor=0):
     N_octaves = int(np.floor(np.log(fs/2) / np.log(2-overlap_factor) - np.log2(16) + 1))
     
     return N_octaves
     
def fractional_octave_filterbank(x, fs, L_DFT, L_block, L_feed, overlap_factor):
    f, t, Pxx = spectrogram(x, fs=fs, window='hann', nperseg=L_block, noverlap=L_block-L_feed, mode='psd')
    
    N_octaves = get_N_bands(L_DFT, fs, fraction=1, overlap_factor=overlap_factor)
    
    retval = np.zeros(shape=(N_octaves, Pxx.shape[1]))
    vec_f_center = np.zeros(N_octaves)
    
    for idx_octave in range(N_octaves):
        idx_lower, idx_center, idx_upper = get_bin_indices_for_octaves(idx_octave, L_DFT, fs, overlap_factor=overlap_factor)
        
        retval[idx_octave,:] = np.sum(Pxx[idx_lower:idx_upper, :], axis=0)
        vec_f_center[idx_octave] = idx_center * fs / L_DFT
        
    return retval, vec_f_center, t
     
    
if __name__ == '__main__':
    # step through all bands and plot the areas

    fs = 44100
    L_DFT = 1024
    
    plt.close('all')
    
    # 1. non-overlapping octave bands
    vec_f = np.arange(L_DFT/2+1) * fs / L_DFT
    
    # calculate the number of octaves for a given set of parameters
    N_octaves = get_N_bands(L_DFT, fs, fraction=1, overlap_factor=0.5)
    
    
    plt.figure(1)
    plt.plot(vec_f, np.zeros(L_DFT/2+1))
    
    ax = plt.gca()
    for idx_octave in range(N_octaves):
        idx_lower, idx_center, idx_upper = get_bin_indices_for_octaves(idx_octave, L_DFT, fs, overlap_factor=0.5)
        
        f_lower = idx_lower * fs/L_DFT
        f_center = idx_center * fs/L_DFT
        f_upper = idx_upper * fs/L_DFT
        
        ax.add_patch(patches.Rectangle([f_lower, 0], f_upper-f_lower, 1-idx_octave/(N_octaves*2), fill=False))
        ax.add_patch(patches.Arrow(f_center, 0, 0, 1))
    
    plt.ylim([-0.1, 1.1])
    plt.xlim((-1000, fs/2+1000))
#    plt.gca().set_xscale("log")
    plt.show()
    
    if False:
        # 2. overlapping octave bands
        overlap_factor = 0.5
        fraction = 1
        
        # calculate the number of octaves for a given set of parameters
        N_octaves = get_N_bands(L_DFT, fs, fraction=fraction, overlap_factor=overlap_factor)
        print(N_octaves)
        
        
        plt.figure(2)
        plt.plot(vec_f, np.zeros(L_DFT/2+1))
        
        ax = plt.gca()
        for idx_octave in range(N_octaves):
            f_lower, f_center, f_upper = get_bin_indices_for_octaves(idx_octave, L_DFT, fs, fraction, overlap_factor)
            
            ax.add_patch(patches.Rectangle([f_lower, 0], f_upper-f_lower, 1-idx_octave/(N_octaves*2), fill=False))
            ax.add_patch(patches.Arrow(f_center, 0, 0, 1))
        
        plt.ylim([-0.1, 1.1])
        plt.xlim((-1000, fs/2+1000))
    #    plt.gca().set_xscale("log")
        plt.show()