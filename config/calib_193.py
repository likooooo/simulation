verbose = -1
wavelength = 0.193 # um

n = 1
tilesize, ambit = [64 * n, 64 * n], [wavelength * 1, wavelength* 1]
maxNA, maxSigma = 1.2, 0
# gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU.gg"
gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_7.ss" #0.23422
# gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_2.ss"
# gauge_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_bug_record.ss" # 0.111771
oas_file = "/home/like/model_data/X_File/LG40_poly_File/LG40_PC_CDU_Contour_Mask_L300.oas"
cell_name, layer_id, dbu = "JDV_M", 300, 0.00025
mask_USF = 8 # -> step=20 dbu (5nm)
mask_edge_dissect_coef = 0.05 # dissect 1dbu
cutline_debug_list = []#list(range(0, 32))

gauss_laguerre_sigma_in_dbu = 30
associate_order = 1
laguerre_order = 2
resist_coefficients, threshold_guess = [], 0.5
# resist_coefficients, threshold_guess = [-0.855563,0.921637,1.87816,-2.86478,2.24506], 0.0533628     

# resist_coefficients, threshold_guess = [-1.33646,0.527619,2.21799,-3.57927,0.0658159], -0.12339 

# resist_coefficients, threshold_guess = [-0.329115,-0.130417,1.56463,-1.20995,-0.362648], 0.231817  
# resist_coefficients, threshold_guess = [-1.33646,0.527619,2.21799,-3.57927,0.0658159]  ,  -0.12339 