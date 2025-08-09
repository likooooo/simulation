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

background = 1
absorber = 0
g2 = geo_manager([[
        [int(-cd/2/grid_info_2d.dbu), int(-cd/2/grid_info_2d.dbu)], 
        [int(-cd/2/grid_info_2d.dbu), int(cd/2/grid_info_2d.dbu)], 
        [int(cd/2/grid_info_2d.dbu), int(cd/2/grid_info_2d.dbu)], 
        [int(cd/2/grid_info_2d.dbu), int(-cd/2/grid_info_2d.dbu)]
    ]
])
mask = binary_mask.create(absorber, background, grid_info_2d, g2.get_vertex())
# test 1
if cd % 2 == 0:
    thin = thin_mask.create(absorber, background, grid_info_2d, g2.get_all_edges())
    error = [abs(a - b) for a, b in zip(thin, mask)]
    assert(max(error) < grid_info_2d.dbu)
# test 2
golden_amp = np.loadtxt(os.path.join(os.path.abspath(os.path.expanduser(golden_root_dir)), "thin_mask_2d_mag_binary.txt"), delimiter='\t').flatten().tolist()
amp_error = [abs(a - abs(x)) for a, x in zip(golden_amp, mask)]
assert(max(amp_error) < grid_info_2d.dbu)
# test 3
nearfield_diffract =  np.fft.fftshift(np.fft.fft2(np.array(mask, dtype=np.complex64).reshape((size, size))))
difract = diffraction(grid_info_2d)
print(difract.update_diffraction_source_points(nearfield_diffract.flatten().tolist()))
pupil_intensity = difract.get_imaging_pupil_intensity([100, 100])
golden_pupil = np.loadtxt(os.path.join(os.path.abspath(os.path.expanduser(golden_root_dir)), "thin_mask_2d_pupil_intensity_binary.txt"), delimiter='\t').flatten().tolist()
amp_error = [a - x for a, x in zip(golden_pupil, pupil_intensity)]
# import matplotlib.pyplot as plt
# plt.imshow(np.array(amp_error).reshape((100, 100)), cmap='viridis')
# plt.show()
assert(max(np.abs(amp_error)) < grid_info_2d.dbu)



