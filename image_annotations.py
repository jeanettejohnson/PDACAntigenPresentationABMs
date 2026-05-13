import pandas as pd
import numpy as np


annotations = pd.read_table("/Users/jeanette.johnson/OneDrive - University of Maryland School of Medicine/JHH387_multipage_ROIs/ROI2Ki67.txt")
annotations[["cell_type"]] = "" 
annotations[["Ki67"]] = False
annotations[["HLA-DR"]] = False


annotations = pd.DataFrame(annotations)

annotations["cell_type"][annotations["Classification"].str.startswith("Tumor")] = "tumor" # ignore the warning, this works
annotations["cell_type"][annotations["Classification"].str.startswith("ductal")] = "ductal"
annotations["cell_type"][annotations["Classification"].str.startswith("CAF")] = "CAF"
annotations["cell_type"][annotations["Classification"].str.startswith("CD4")] = "CD4 T cell"
annotations["cell_type"][annotations["Classification"].str.startswith("CD8")] = "CD8 T cell"
annotations["cell_type"][annotations["Classification"].str.startswith("Other")] = "other_cell"



annotations["HLA-DR"][annotations["Classification"].str.contains("HLA-DR")] = True
annotations["Ki67"][annotations["Classification"].str.contains("Ki67")] = True

annotations.to_csv("/Users/jeanette.johnson/OneDrive - University of Maryland School of Medicine/JHH387_multipage_ROIs/ROI2_Ki67_HLADR_annotations.csv")