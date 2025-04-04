verbose = 0
wavelength = 0.193 # um

tilesize, ambit = [64, 64], [wavelength * 1, wavelength* 1]
maxNA, maxSigma = 1.2, 0
gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU.gg"
# gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_2.ss"
# gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_bug_record.ss"
oas_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_Contour_Mask_L300.oas"
cell_name, layer_id, dbu = "JDV_M", 300, 0.00025
mask_USF = 8 # -> step=20 dbu (5nm)
mask_edge_dissect_coef = 0.05 # dissect 1dbu
threshold_guess =  0.264401
cutline_debug_list =[]# list(range(0, 32))