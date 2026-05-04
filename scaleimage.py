import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler


image = pd.read_csv("/Users/jeanette.johnson/Downloads/ROI1Collagen.csv")


scaler = MinMaxScaler()

image[["ecm"]] = scaler.fit_transform(image[["ecm"]])*10

# squish and then 

x = image["x"].values
y = image["y"].values
if len(x) > 1: # I know, right, what ST data set will have one spot??
    xL = np.min(x)
    xR = np.max(x)
    yL = np.min(y)
    yR = np.max(y)
else:
    print("Single ST spot?")
    xL, xR = x[0] + [-0.5,0.5]
    yL, yR = y[0] + [-0.5,0.5]

spatial_factors = [1600 / (xR-xL), 1600 / (yR-yL)] # factors for scaling each dimension to an interval of length 1


spatial_factor = min(spatial_factors)
width = (xR-xL)*spatial_factor
height = (yR-yL)*spatial_factor

x0 = 0.5*(-800+800 - width)
y0 = 0.5*(-800+800 - height)

spatial_base_coords = image[["x", "y"]] - [xL,yL]

spatial_base_coords = spatial_base_coords * [spatial_factor,spatial_factor]
spatial_base_coords = spatial_base_coords - [800, 800]
spatial_base_coords[["z"]] = 0
spatial_base_coords[["ecm"]] = image[["ecm"]]

# spatial_base_coords[["y"]] = 1- spatial_base_coords[["y"]]
# y = 1-y
# flip image twice so it matches with BIWT 
coordarray = np.array(spatial_base_coords)
#coordarray = np.flip(coordarray, axis=1)
#coordarray = np.flip(coordarray, axis=0)
#coordarray = np.rot90(coordarray, k=2)
spatial_base_coords = pd.DataFrame(np.array(coordarray))

spatial_base_coords.to_csv("/Users/jeanette.johnson/OneDrive - University of Maryland School of Medicine/JHH387_multipage_ROIs/ROI001_ROI_001/ROI1Collagen.csv")