import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathDetectionObject

// ── KNOBS ──────────────────────────────────────────────────────────────────
def separator = "\t"
// ───────────────────────────────────────────────────────────────────────────

def project   = getProject()
def exportDir = new File("/Users/jeanette.johnson/PDACAntigenPresentationABMs/qupath_detections")
exportDir.mkdirs()

def entries = project.getImageList()
print "Found ${entries.size()} images"

entries.each { entry ->
    def imageName  = entry.getImageName()
    def outputFile = new File(exportDir, imageName + "_detections.txt")

    def imageData   = entry.readImageData()
    def nDetections = imageData.getHierarchy().getDetectionObjects().size()

    if (nDetections == 0) {
        print "Skipping ${imageName} — no detections"
        return
    }

    new MeasurementExporter()
        .imageList([entry])
        .separator(separator)
        .exportType(PathDetectionObject.class)
        .exportMeasurements(outputFile)

    print "Exported ${nDetections} detections: ${imageName}"
}

print "Done — files saved to: ${exportDir}"
