import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathDetectionObject

// ── KNOBS ──────────────────────────────────────────────────────────────────
def separator = "\t"
// ───────────────────────────────────────────────────────────────────────────

def exportDir = new File("/Users/jeanette.johnson/PDACAntigenPresentationABMs/qupath_detections")
exportDir.mkdirs()

def entry      = getProjectEntry()
def imageName  = entry.getImageName()
def outputFile = new File(exportDir, imageName + "_detections.txt")

def imageData   = getCurrentImageData()
def nDetections = imageData.getHierarchy().getDetectionObjects().size()

if (nDetections == 0) {
    print "No detections found in ${imageName} — nothing exported"
    return
}

// Flush in-memory state (including latest classifier results) to disk
entry.saveImageData(imageData)
print "Saved image data for ${imageName}"

new MeasurementExporter()
    .imageList([entry])
    .separator(separator)
    .exportType(PathDetectionObject.class)
    .exportMeasurements(outputFile)

print "Exported ${nDetections} detections to: ${outputFile}"
