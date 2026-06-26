import java.io.BufferedReader;
import java.io.FileReader;
import qupath.lib.objects.PathAnnotationObject;
import qupath.lib.roi.RectangleROI;

def imageData = getCurrentImageData();

// Create BufferedReader
// Source - https://stackoverflow.com/a/9084256
// Posted by Michael Boselowitz


String[] channels = new String[56];
int i = 0;
try (BufferedReader br = new BufferedReader(new FileReader("/Users/jeanette.johnson/OneDrive - University of Maryland School of Medicine/JHH_IMC/multipage/JHH387_multipage_ROIs/ROI001_ROI_001/ROI001_ROI_001_summary.txt"));
) {
    String line;
    br.readLine()
    while ((line = br.readLine()) != null) {
        String[] data = line.split("\t");
        print(data[2]);
        print(line)
        channels[i] = data[2];
        i=i+1;
    }
} catch (IOException e) {
    e.printStackTrace();
}

print(channels)

setChannelNames(channels)