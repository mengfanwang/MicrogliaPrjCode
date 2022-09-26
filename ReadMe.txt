Files, functions, inputs and outputs.

ReadFile: Convert czi file to mat and visiualization.
.czi -> .mat data.tif

VarianceStabilization: stabilize the variances.
	histogramCount: count histogram
	convexOptimization: find transform function
data -> stabilizeFunction variance stabilization.tif

ForegroundDetection: detect foreground
	growRegion
data stabilizeFunction variance -> foreground

SomaDetection: detect soma
foreground -> soma

Registration: registrate data based on soma
foreground soma -> registration

CellSepration: cut one cell and build skeleton
	treeBuilding
foreground soma -> cell.tif result graph(Node/G) cell(isolated region) cell_bg(with background)

GraphMatching: match adjacent frames by min-cost flow
graph -> matchingResult

Matchibg2Connection: find over/under connection
matchingResult graph result cell cell_bg -> result_post