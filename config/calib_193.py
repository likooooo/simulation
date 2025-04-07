verbose = -1
wavelength = 0.193 # um

n = 1
tilesize, ambit = [64 * n, 64 * n], [wavelength * 1, wavelength* 1]
maxNA, maxSigma = 1.2, 0
gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU.gg"
# gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_2.ss" #0.23422
# gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_bug_record.ss" # 0.111771
oas_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_Contour_Mask_L300.oas"
cell_name, layer_id, dbu = "JDV_M", 300, 0.00025
mask_USF = 8 # -> step=20 dbu (5nm)
mask_edge_dissect_coef = 0.05 # dissect 1dbu
threshold_guess = 0.5
cutline_debug_list = [0]#list(range(0, 32))