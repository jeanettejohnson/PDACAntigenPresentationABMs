from imctools.io.mcd.mcdparser import McdParser

fn_mcd = "/Users/jeanette.johnson/OneDrive - University of Maryland School of Medicine/JHH_IMC/JHH417R/JHH417R/JHH417R.mcd"

parser = McdParser(fn_mcd)

# Get original metadata in XML format
xml = parser.get_mcd_xml()

# Get parsed session metadata (i.e. session -> slides -> acquisitions -> channels, panoramas data)
session = parser.session

# Get all acquisition IDs
ids = parser.session.acquisition_ids

# The common class to represent a single IMC acquisition is AcquisitionData class.
# Get acquisition data for acquisition with id 2
ac_data = parser.get_acquisition_data(3)


# imc acquisitions can yield the image data by name (tag), label or index
channel_image1 = ac_data.get_image_by_name('Ir191')
channel_image2 = ac_data.get_image_by_label('Histone_phospho_125((2468))Eu153')
channel_image3 = ac_data.get_image_by_index(7)

# or can be used to save OME-TIFF files
fn_out ='/Users/jeanette.johnson/OneDrive - University of Maryland School of Medicine/JHH_IMC/JHH417R.ome.tiff'
ac_data.save_ome_tiff(fn_out, names=['Ir191', 'Yb172'])


from imctools.converters import mcdfolder_to_imcfolder

mcdfolder_to_imcfolder("/Users/jeanette.johnson/JHHPDACIMCAnalysis", "/Users/jeanette.johnson/JHHPDACIMCAnalysis/ometiff")
# save multiple standard TIFF files in a folder
ac_data.save_tiffs("/home/anton/tiffs", compression=0, bigtiff=False)

# as the mcd object is using lazy loading memory maps, it needs to be closed
# or used with a context manager.
parser.close()



