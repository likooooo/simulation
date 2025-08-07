#!/bin/bash
# same behavior check
golden_output="
# create grid info in default mode.
* grid info                    
---------------------------------
tilesize      :  [256,1]       
spatial start :  [0,-200]      
spatial step  :  [1.5625,0]    
fourier start :  [0,0]         
fourier step  :  [0.4825,inf]  
coords place  :  |*--|corner
# create grid info in bloch mode.
* grid info                    
---------------------------------
tilesize      :  [256,1]       
spatial start :  [0,-200]      
spatial step  :  [1.5625,0]    
fourier start :  [0,0]         
fourier step  :  [0.4825,inf]  
coords place  :  |*--|corner
# create grid info in OPC mode.
* grid info                           
----------------------------------------
tilesize      :  [256,1]              
spatial start :  [-4946.67,-220.104]  
spatial step  :  [40.2083,40.2083]    
fourier start :  [0,0]                
fourier step  :  [0.01875,4.8]
coords place  :  |*--|corner
"
for (( i=0; i<=2; i++ )); do
    output_cmd=$(./test_grid_info $i  256 1 193 0 1.2  0 -200 400 -200)
    if ! echo "$output_cmd" | grep -Fq "$golden_output"; then
        echo "error output"
        echo "$output_cmd"
        exit 1
    fi
    diff <(echo "$output_cmd") <(python3 -c "from grid_info import *;a = create_grid_info($i, [256, 1], 193, 0, 1.2, [[0, -200], [400,-200]], 1e-6)") || exit $?
done

# Wafer P.O.V
# fourier step : 0.902
# cutoff freq  : 1
# max freq     : 1/NA
# 目前 sigma 没有对应 hyper litho 的功能，建议设成 0
./test_grid_info 1  1 1 13  0 0.9  -8 0 8 0

# Mask P.O.V
# mask POV 的 tilesize, 和 fourier step 的变化，等效于在 Wafer P.O.V 空域边界填 0，填充至 Demangnification 倍大小