import pandas as pd

qupath_anno = pd.read_table("/Users/jeanette.johnson/OneDrive - University of Maryland School of Medicine/JHH_IMC/annotation_coordinates/JHH387ROI3.txt")
qupath_anno = pd.DataFrame(qupath_anno)

qupath_for_physicell = qupath_anno[["Classification","Centroid X µm", "Centroid Y µm"]]

qupath_for_physicell.columns = ["type", "x", "y"]

qupath_for_physicell.to_csv("/Users/jeanette.johnson/OneDrive - University of Maryland School of Medicine/JHH_IMC/annotation_coordinates/JHH387ROI3_for_physicell.csv", index=False)
