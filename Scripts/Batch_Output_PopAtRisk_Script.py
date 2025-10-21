import csv
from qgis.core import QgsProject, QgsVectorLayer, QgsWkbTypes

# Prepare output CSV path (update if needed)
output_csv = "/Users/la-lcouper/Downloads/zstat_sums6.csv"  # <-- Update path if needed

# Initialize list to hold rows
output_data = [("LayerName", "ZStat_Sum")]

# Loop through layers
for layer in QgsProject.instance().mapLayers().values():
    if not isinstance(layer, QgsVectorLayer):
        continue
    if layer.geometryType() != QgsWkbTypes.PolygonGeometry:
        continue
    if not layer.fields().indexFromName("zstat_sum") >= 0:
        print(f"⚠️  No 'zstat_sum' field in layer: {layer.name()}")
        continue

    # Sum all zstat_sum values across features
    total_sum = 0
    for feature in layer.getFeatures():
        val = feature["zstat_sum"]
        if val is not None:
            total_sum += val

    output_data.append((layer.name(), total_sum))

# Write to CSV
with open(output_csv, mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerows(output_data)

print(f"✅ Exported {len(output_data)-1} layers to: {output_csv}")