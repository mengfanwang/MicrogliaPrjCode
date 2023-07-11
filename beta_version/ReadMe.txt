Files, functions, inputs and outputs.

ReadFile: Convert czi file to mat and visiualization.
.czi -> .mat data.tif

Registration: registrate vessel data
data1 -> shifts1
vessel_OrderStatisticThresholding_v6: detect vessel and ablation site

VarianceStabilization: stabilize the variances.
	histogramCount: count histogram
	convexOptimization: find transform function
data -> stabilizeFunction variance stabilization.tif

vessel_detection: vessel_OrderStatisticThresholding_v6
registration -> foreground_c1
gaussian_sigma is a important variable

glia ForegroundDetection: detect foreground
	growRegion
data stabilizeFunction variance -> foreground

soma detection

distance_probability:
foreground_c2 -> distance_map 

distance_maximum

output