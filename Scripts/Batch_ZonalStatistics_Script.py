from qgis.analysis import QgsZonalStatistics
from qgis.core import QgsProject, QgsVectorLayer, QgsWkbTypes

# Get the loaded raster layer by name
raster_layer = QgsProject.instance().mapLayersByName("CAPOP_2020_100m_TOTAL")[0]

# Loop through all layers in the project
for layer in QgsProject.instance().mapLayers().values():
    if not isinstance(layer, QgsVectorLayer):
        continue
    if layer.geometryType() != QgsWkbTypes.PolygonGeometry:
        continue

    print(f"Running zonal statistics on layer: {layer.name()}")

    # Run zonal statistics (using positional arguments)
    zs = QgsZonalStatistics(
        layer,                  # vector layer
        raster_layer,           # raster layer
        "zstat_",               # prefix for new fields
        1,                      # raster band (usually 1)
        QgsZonalStatistics.Sum  # statistic to calculate
    )
    zs.calculateStatistics(None)

print("✅ Zonal statistics complete.")
