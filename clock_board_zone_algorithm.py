from qgis.core import (
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterRasterLayer,
    QgsFeature,
    QgsGeometry,
    QgsPointXY,
    QgsFields,
    QgsField,
    QgsWkbTypes,
    QgsProcessingException,
    QgsProcessing,
    QgsFeatureSink,
    QgsVectorLayer,
    QgsCoordinateTransform,
    QgsCoordinateReferenceSystem
)
from qgis.PyQt.QtCore import QVariant, QCoreApplication
import math

class ClockBoardZoneAlgorithm(QgsProcessingAlgorithm):
    def initAlgorithm(self, config=None):
        # Input point layer (multiple centers)
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                'CENTER_LAYER',
                'Center Points Layer (point layer; all features used)',
                [QgsProcessing.TypeVectorPoint]
            )
        )
        # Optional raster for zonal statistics
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                'RASTER',
                'Input Raster Layer for Zonal Statistics',
                optional=True
            )
        )
        # Number of circles
        self.addParameter(
            QgsProcessingParameterNumber(
                'CIRCLE_COUNT',
                'Number of Circles',
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=3
            )
        )
        # Circle size in meters
        self.addParameter(
            QgsProcessingParameterNumber(
                'CIRCLE_SIZE',
                'Circle Size (meters)',
                type=QgsProcessingParameterNumber.Double,
                defaultValue=1000
            )
        )
        # Segments per circle
        self.addParameter(
            QgsProcessingParameterNumber(
                'SEGMENT_COUNT',
                'Number of Segments',
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=8
            )
        )
        # Output vector
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                'OUTPUT',
                'Output Zones with Zonal Statistics',
                type=QgsProcessing.TypeVectorPolygon
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        # Read all center points
        source = self.parameterAsSource(parameters, 'CENTER_LAYER', context)
        centers = list(source.getFeatures())
        if not centers:
            raise QgsProcessingException('No center points provided.')

        circle_count = self.parameterAsInt(parameters, 'CIRCLE_COUNT', context)
        circle_size = self.parameterAsDouble(parameters, 'CIRCLE_SIZE', context)
        segment_count = self.parameterAsInt(parameters, 'SEGMENT_COUNT', context)

        # Prepare fields: center_id, circle, segment, zone_id
        fields = QgsFields()
        fields.append(QgsField('center_id', QVariant.Int))
        fields.append(QgsField('circle', QVariant.Int))
        fields.append(QgsField('segment', QVariant.Int))
        fields.append(QgsField('zone_id', QVariant.String))

        # Create memory layer in EPSG:3857
        target_crs = QgsCoordinateReferenceSystem('EPSG:3857')
        mem_layer = QgsVectorLayer(f"Polygon?crs={target_crs.authid()}", "zones", "memory")
        prov = mem_layer.dataProvider()
        prov.addAttributes(fields.toList())
        mem_layer.updateFields()

        # Helper: compute number of arc points per radius
        def calc_arc_points(r):
            pts = int((2 * math.pi * r) / 50.0)
            return max(8, min(pts, 200))

        all_feats = []
        # Loop through centers
        for cid, feat in enumerate(centers, start=1):
            geom = feat.geometry()
            crs_src = source.sourceCrs()
            xform = QgsCoordinateTransform(crs_src, target_crs, context.transformContext())
            geom.transform(xform)
            center_pt = geom.asPoint()

            # Build rings for this center
            for i in range(circle_count):
                inner_r = i * circle_size
                outer_r = (i + 1) * circle_size
                arcs = calc_arc_points(outer_r)

                for j in range(segment_count):
                    ang0 = (j * 2 * math.pi) / segment_count
                    ang1 = ((j + 1) * 2 * math.pi) / segment_count
                    pts = []
                    # Outer arc
                    for k in range(arcs + 1):
                        a = ang0 + (k / arcs) * (ang1 - ang0)
                        pts.append(QgsPointXY(
                            center_pt.x() + outer_r * math.cos(a),
                            center_pt.y() + outer_r * math.sin(a)
                        ))
                    # Inner arc
                    if inner_r == 0:
                        pts.append(center_pt)
                    else:
                        for k in range(arcs + 1):
                            a = ang1 - (k / arcs) * (ang1 - ang0)
                            pts.append(QgsPointXY(
                                center_pt.x() + inner_r * math.cos(a),
                                center_pt.y() + inner_r * math.sin(a)
                            ))
                    # Close polygon
                    if pts[0] != pts[-1]:
                        pts.append(pts[0])

                    f = QgsFeature()
                    f.setFields(fields)
                    zid = f"{cid}_{i+1}.{j+1}"
                    f.setAttributes([cid, i+1, j+1, zid])
                    f.setGeometry(QgsGeometry.fromPolygonXY([pts]))
                    all_feats.append(f)

                    if feedback.isCanceled():
                        break
                if feedback.isCanceled():
                    break
            if feedback.isCanceled():
                break

        # Dissolve circle=1 per center
        final = []
        for cid in range(1, len(centers)+1):
            # circle 1 features for this center
            c1 = [f.geometry() for f in all_feats if f['center_id']==cid and f['circle']==1]
            if c1:
                dg = QgsGeometry.unaryUnion(c1)
                nf = QgsFeature()
                nf.setFields(fields)
                nf.setAttributes([cid, 1, 1, f"{cid}_1"])
                nf.setGeometry(dg)
                final.append(nf)
            # add non-circle-1
            final.extend([f for f in all_feats if not (f['center_id']==cid and f['circle']==1)])

        prov.addFeatures(final)
        mem_layer.updateExtents()

        # Zonal stats if raster provided
        raster = self.parameterAsRasterLayer(parameters, 'RASTER', context)
        if raster:
            from qgis.analysis import QgsZonalStatistics
            zs = QgsZonalStatistics(mem_layer, raster, 'zone_', 1,
                                    QgsZonalStatistics.Sum | QgsZonalStatistics.Mean)
            zs.calculateStatistics(None)

        # Prepare sink
        (sink, sid) = self.parameterAsSink(
            parameters, 'OUTPUT', context,
            mem_layer.fields(), mem_layer.wkbType(), target_crs
        )
        for f in mem_layer.getFeatures():
            sink.addFeature(f, QgsFeatureSink.FastInsert)

        return {'OUTPUT': sid}

    def name(self):
        return 'clock_board_zone'

    def displayName(self):
        return self.tr('Clock Board Zone (Multi-Center)')

    def group(self):
        return self.tr('Urban and Regional Planning Analysis')

    def groupId(self):
        return 'urbanremotesensing'

    def createInstance(self):
        return ClockBoardZoneAlgorithm()

    def tr(self, s):
        return QCoreApplication.translate('ClockBoardZoneAlgo', s)

    def shortHelpString(self):
        return """
        <b>Clock Board Zone</b><br>
        This tool generates circular buffer zones from input points for spatial performance analysis.Users can specify the number of circles and their radii. For the first circle, all segments are dissolved into a single polygon before zonal statistics (sum and mean) are calculated based on an input raster layer (optional).<br><br>
        <b>Purpose:</b><br>
        This tool is designed to support spatial performance analysis by creating multiple buffer zones around key locations and overlaying them with raster data for extracting relevant statistics.<br><br>
        <b>Credits:</b><br>
        This tool was developed by <b>Firman Afrianto</b> as part of the <b>Digital Data Processing</b> course in the <b>Interdisciplinary Master's Program in Urban and Regional Planning, Faculty of Engineering, University of Indonesia (FTUI)</b>. 
        This course is taught by:
        <ul>
            <li>Hendricus Andy Simarmata</li>
            <li>Firman Afrianto</li>
            <li>Jati Pratomo</li>
        </ul>
        <br>
        <b>Note:</b> It is recommended to use QGIS 3.38 and CRS Pseudomercator (EPSG:3857) for optimal performance.
        """
