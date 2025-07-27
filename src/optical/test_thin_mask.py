#!/usr/bin/env python3
from grid_info import *
from math import *
from cmath import exp, phase
import numpy as np
import os

def polar_to_complex(amplitude, phase):
    return amplitude * exp(1j * phase)

golden_root_dir = "~/doc/github/simulation/resource/golden_data_latest/near_field/hyper_lith_kirchhoff/"
wavelength = 13.0
NA = 0.9
pitch = 16
cd = 8
dx = 1
size = int(ceil(pitch / dx))

grid_info_2d = grid_info_2d_s.create_grid_info_bloch_mode([size, size], wavelength, 0.0, NA, [[-pitch/2, -pitch/2], [pitch/2, pitch/2]], 1e-6)
print(grid_info_2d)

background = polar_to_complex(0.86065, 0)
absorber = polar_to_complex(0.0968, -2.675)
g2 = geo_manager([[
        [int(-cd/2/grid_info_2d.dbu), int(-cd/2/grid_info_2d.dbu)], 
        [int(-cd/2/grid_info_2d.dbu), int(cd/2/grid_info_2d.dbu)], 
        [int(cd/2/grid_info_2d.dbu), int(cd/2/grid_info_2d.dbu)], 
        [int(cd/2/grid_info_2d.dbu), int(-cd/2/grid_info_2d.dbu)]
    ]
])
mask = binary_mask.create(absorber, background, grid_info_2d, g2.get_vertex())
if cd % 2 == 0:
    thin = thin_mask.create(absorber, background, grid_info_2d, g2.get_all_edges())
    error = [abs(a - b) for a, b in zip(thin, mask)]
    assert(max(error) < grid_info_2d.dbu)

golden_amp = np.loadtxt(os.path.join(os.path.abspath(os.path.expanduser(golden_root_dir)), "thin_mask_2d_mag.txt"), delimiter='\t').flatten().tolist()
amp_error = [abs(a - abs(x)) for a, x in zip(golden_amp, mask)]
assert(max(amp_error) < grid_info_2d.dbu)

golden_phase = np.loadtxt(os.path.join(os.path.abspath(os.path.expanduser(golden_root_dir)), "thin_mask_2d_phase.txt"), delimiter='\t').flatten().tolist()
phase_error = [abs(a - phase(x)) for a, x in zip(golden_phase, mask)]
assert(max(phase_error) < grid_info_2d.dbu)



