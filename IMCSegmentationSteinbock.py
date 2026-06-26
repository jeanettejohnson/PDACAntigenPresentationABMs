import steinbock as steinbock
steinbock.preprocessing.preprocess_image("/Users/jeanette.johnson/HTANWUSTL/IMC/HT168P1_S1H3_L1_L4_ROI_002.ome.tiff")


steinbock.create_image_info("/Users/jeanette.johnson/HTANWUSTL/IMC/HT168P1_S1H3_L1_L4_ROI_002.ome.tiff")

# Import Python libraries
# This code uses python v3.12.0, tifffile v2023.9.26, matplotlib v3.8.2
import tifffile
import matplotlib.pyplot as plt

# Option 1: Load full resolution image channels
# The following may produce a warning: 'OME series cannot read multi-file pyramids'. This is because tifffile does not support loading a pyramidal multi-file OME-TIFF file. Only the full resolution (level=0) data will load for all channels in the directory.
fullres_multich_img = tifffile.imread(
"/Users/jeanette.johnson/HTANWUSTL/IMC/data/ometiff/HT168P1_S1H3_L1_L4_ROI_002.ome.tiff", is_ome=True, level=0, aszarr=False)

# Examine shape of array (number of channels, height, width), e.g. (4, 40867, 31318)
fullres_multich_img.shape

# Extract number of channels, e.g. 4
n_ch = fullres_multich_img.shape[0]

# Plot each channel
fig, axes = plt.subplots(ncols=n_ch, nrows=1, squeeze=False, figsize=plt.figaspect(1/2))
for i in range(n_ch):
    axes[0, i].imshow(fullres_multich_img[i], cmap="gray")
    axes[0, i].set_title(f"Channel: {i}")
plt.savefig('tifffile_fullres_four_channels.png')

# Option 2: Load a single channel image at any resolution, e.g., level=0 or level=1. Note 'is_ome' is set to False.
# Load one of the multi-file OME-TIFF files as a regular TIFF file at full resolution.
fullres_img_tiff = tifffile.imread(
    "/Users/jeanette.johnson/HTANWUSTL/IMC/data/ometiff/HT168P1_S1H3_L1_L4_ROI_002.ome.tiff", is_ome=False, level=0)

# Now load the file at downsampled resolution
downsampled_img = tifffile.imread(
    "/Users/jeanette.johnson/HTANWUSTL/IMC/data/ometiff/HT168P1_S1H3_L1_L4_ROI_002.ome.tiff", is_ome=False, level=1)
# Examine shape of array (number of channels, height, width)
downsampled_img.shape

# Plot the full resolution and downsampled images side-by-side
fig, axes = plt.subplots(ncols=2, nrows=1, squeeze=False)
axes[0, 0].imshow(fullres_img_tiff, cmap="gray")
axes[0, 0].set_title(f"Full resolution: {fullres_img_tiff.shape}")
axes[0, 1].imshow(downsampled_img, cmap="gray")
axes[0, 1].set_title(f"Downsampled: {downsampled_img.shape}")
plt.savefig('example_fullres_downsample.png')
