'''
Copyright 2021 Esri

Licensed under the Apache License Version 2.0 (the "License");
you may not use this file except in compliance with the License.

You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

See the License for the specific language governing permissions and
limitations under the License.
'''

# -*- coding: utf-8 -*-
import arcpy, os
import math
import numpy as np
from typing import NamedTuple

##******************************************##
##               Main Class                 ##
##******************************************##

class IdentifySparseVertices(object):
    def __init__(self):
        self.label = "Identify Sparse Vertices"
        self.description = "This tool takes an input polyline feature class that is configured with GCS and calculates segments which have vertices that are too far apart. The output feature class will contain a record for each segment that exceeds that maximum distance based on the current XY/Z tolerances."
        self.canRunInBackground = False

    def getParameterInfo(self):

        params = []

        # Input Feature Class.
        param0 = arcpy.Parameter(
            displayName="Input Feature Class",
            name="in_feature",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input")
        param0.filter.list = ["Polyline"]

        # Output Feature Class.
        param1 = arcpy.Parameter(
            displayName="Output Feature Class",
            name="out_feature",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Output")

        params = [param0, param1]

        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        
        # Get Parameter info.
        in_feature = parameters[0].valueAsText
        out_feature = parameters[1].valueAsText

        # Clear previous errors.
        parameters[0].clearMessage()
        parameters[1].clearMessage()

        # Validate input feature class.
        if parameters[0].altered:

            result = ValidateFeatureClass(in_feature)
            if result:
                parameters[0].setErrorMessage(result)
                return

        # Validate output path.
        if parameters[1].altered:

            result = ValidateOutputFeatureClass(out_feature)
            if result:
                parameters[1].setErrorMessage(result)
                return

        return

    def execute(self, parameters, messages):
        
        # Input Feature Class.
        in_feature_class = parameters[0].valueAsText
        sr = arcpy.Describe(in_feature_class).spatialReference

        # Output workspace
        out_workspace = parameters[1].valueAsText

        # Get max distances for vairous tolerances.
        max_distances = GetMaxDistances(in_feature_class, sr)

        # Create the output feature class.
        out_feature_class = CreateOutputFeatureClass(parameters[1].valueAsText, in_feature_class, sr)
        
        # Query sparse vertices.
        record_log = GetSparseVertices(in_feature_class, out_feature_class, max_distances, sr)

        # Write to log file.
        WriteLogFile(record_log, max_distances, sr, out_workspace)

        return

##******************************************##
##                 Methods                  ##
##******************************************##

def GetSparseVertices(in_feature_class, out_feature_class, max_distances, sr):
    
    # Get feature count for progessor.
    featureCount = int(arcpy.GetCount_management(in_feature_class)[0])
    arcpy.SetProgressor("step", "Checking vertex spacing...", 0, featureCount, 1)

    # Set counters to 0.
    XYCount = ZCount = ftCount = mm1Count = cm1Count = cm2Count = cm3Count = cm10Count = 0

    # Get only editable fields. Exclude geometry as we will use a tag to get the full geometry properties.
    fields = arcpy.ListFields(in_feature_class)
    fieldNames = []
    for field in fields:
        if field.type == 'Geometry' or field.editable == False:
            continue
        fieldNames.append(field.name)

    # Create the insert cursor. We will insert each segment that spans more than the max distance.
    with arcpy.da.InsertCursor(out_feature_class, ['Shape@', fieldNames]) as insertCursor:

        # Create the search cursor using the same fields from the insert cursor.
        with arcpy.da.SearchCursor(in_feature_class, [insertCursor.fields]) as searchCursor:
            
            for row in searchCursor:
                
                # Reset flags.
                exceedsXY = exceedsZ = exceedsFt = exceeds01mm = exceeds1mm = exceeds1cm = exceeds2cm = exceeds3cm = exceeds10cm = False

                # Shapes to insert for this row.
                shapes = []

                # Check for null geometry.
                if row[0] is None:
                    arcpy.SetProgressorPosition()
                    continue

                # Iterate each part.
                for part in row[0]:
            
                    # Rest after each part.
                    prevPointGeom = None
                    prevPoint = None

                    # Iterate each vertex within part.
                    for currPoint in part:
                
                        # Use PointGeometry to utilize the distanceTo method.
                        currPointGeom = arcpy.PointGeometry(currPoint, sr)

                        # First point in the polyline
                        if prevPointGeom is None:
                            prevPointGeom = currPointGeom
                            prevPoint = currPoint
                            continue

                        # Reset flag for each pass. Each shape will only consist of 2 points.
                        insert = False

                        # Get distance. Ignore angle.
                        angle,distance = prevPointGeom.angleAndDistanceTo(currPointGeom)

                        # Check if we need to write to output feature class.
                        if distance > max_distances.PrimaryMaxDistance:
                            insert = True

                        # Check if we need to log which max distance these vertices exceed.
                        if distance > max_distances.XYMaxDistance:
                            exceedsXY = True
                        if distance > max_distances.ZMaxDistance:
                            exceedsZ = True
                        if distance > max_distances.FtMaxDistance:
                            exceedsFt = True
                        if distance > max_distances.Mm1MaxDistance:
                            exceeds1mm = True
                        if distance > max_distances.Cm1MaxDistance:
                            exceeds1cm = True
                        if distance > max_distances.Cm2MaxDistance:
                            exceeds2cm = True
                        if distance > max_distances.Cm3MaxDistance:
                            exceeds3cm = True
                        if distance > max_distances.Cm10MaxDistance:
                            exceeds10cm = True

                        # Create shape. We will do the insert after we finished processing the row.
                        if insert:
                            array = arcpy.Array()
                            array.add(prevPoint)
                            array.add(currPoint)
                            shapes.append(arcpy.Polyline(array, sr, True , True ))

                        # Update previous points.
                        prevPointGeom = currPointGeom
                        prevPoint = currPoint

                # Write each shape to the output feature class. 
                for segment in shapes:
                    rowList = list(row)
                    rowList[0] = segment
                    insertCursor.insertRow(tuple(rowList))

                # Update log counters. We only want to increment the counters once per row.
                if exceedsXY:
                    XYCount += 1
                if exceedsZ:
                    ZCount += 1
                if exceedsFt:
                    ftCount += 1
                if exceeds1mm:
                    mm1Count += 1
                if exceeds1cm:
                    cm1Count += 1
                if exceeds2cm:
                    cm2Count += 1                
                if exceeds3cm:
                    cm3Count += 1                
                if exceeds10cm:
                    cm10Count += 1

                # Update progressor.
                arcpy.SetProgressorPosition()

    return RecordLog(featureCount, XYCount, ZCount, ftCount, mm1Count, cm1Count, cm2Count, cm3Count, cm10Count)

def WriteLogFile(record_log, max_distances, sr, workspace):
    
    # Set progressor label.
    arcpy.SetProgressorLabel("Writing to output log file...")

    # Create output at the same location as the db.
    outputDir = os.path.dirname(workspace)
    outputDir = outputDir.rsplit('\\', 1)[0]
    outputLogName = 'SparseVerticesOutput' + r".log"
    outputFileLog = os.path.join(outputDir, outputLogName)
    txtFile = open(outputFileLog, 'w')

    # Write feature count.
    txtFile.write('Input feature class records: %s\n' % record_log.FeatureCount);
    txtFile.write('\n')

    if (record_log.FeatureCount > 0):
        xyInMeters = sr.XYTolerance * (20015077 / 180 )
        txtFile.write('Current XY-Tolerance: %s degree\n' % np.format_float_positional(sr.XYTolerance, trim='-'))
        txtFile.write('Current XY-Tolerance in meters: %s m\n' % xyInMeters)
        txtFile.write('Current Z-Tolerance: %s m\n' % sr.ZTolerance)
        txtFile.write('\n')

        if max_distances.PrimaryMaxDistance == max_distances.XYMaxDistance:
            # Write For XY-tolerance.
            txtFile.write('Based on current tolerance:\n')
            txtFile.write('Max distance in meters: %s m\n' % round(max_distances.XYMaxDistance, 2))
            txtFile.write('Records that contain segments that exceed max distance: %s\n' % record_log.ExceedsXYCount)
            txtFile.write('Percentage of input records that exceed max distance: %s%%\n' % round(((record_log.ExceedsXYCount / record_log.FeatureCount) * 100 ), 2))
            txtFile.write('\n\n')

        else:
            # Write for Z-tolerance.
            txtFile.write('Based on current tolerance:\n')
            txtFile.write('Max Distance: %s m\n' % round(max_distances.ZMaxDistance, 2))
            txtFile.write('Records that contain segments that exceed max distance: %s\n' % record_log.ExceedsZCount)
            txtFile.write('Percentage of input records that exceed max distance: %s%%\n' % round(((record_log.ExceedsZCount / record_log.FeatureCount) * 100 ), 2))
            txtFile.write('\n\n')

        if max_distances.PrimaryMaxDistance != max_distances.FtMaxDistance:
            # Write for tolerance of 0.01 feet.
            txtFile.write('Based on 0.001 ft tolerance \n')
            txtFile.write('Max Distance: %s m\n' % round(max_distances.FtMaxDistance, 2))
            txtFile.write('Records that contain segments that exceed max distance: %s\n' % record_log.ExceedsFtCount)
            txtFile.write('Percentage of input records that exceed max distance: %s%%\n' % round(((record_log.ExceedsFtCount / record_log.FeatureCount) * 100 ), 2))
            txtFile.write('\n')

        if max_distances.PrimaryMaxDistance != max_distances.Mm1MaxDistance:
            # Write for tolerance of 0.001 meters.
            txtFile.write('Based on 0.001 m tolerance \n')
            txtFile.write('Max Distance: %s m\n' % round(max_distances.Mm1MaxDistance, 2))
            txtFile.write('Records that contain segments that exceed max distance: %s\n' % record_log.Exceeds1MmCount)
            txtFile.write('Percentage of input records that exceed max distance: %s%%\n' % round(((record_log.Exceeds1MmCount / record_log.FeatureCount) * 100 ), 2))
            txtFile.write('\n')

        if max_distances.PrimaryMaxDistance != max_distances.Cm1MaxDistance:
            # Write for tolerance of 1 cm.
            txtFile.write('Based on 0.01 m tolerance\n')
            txtFile.write('Max Distance: %s m\n' % round(max_distances.Cm1MaxDistance, 2))
            txtFile.write('Records that contain segments that exceed max distance: %s\n' % record_log.Exceeds1CmCount)
            txtFile.write('Percentage of input records that exceed max distance: %s%%\n' % round(((record_log.Exceeds1CmCount / record_log.FeatureCount) * 100 ), 2))
            txtFile.write('\n')

        if max_distances.PrimaryMaxDistance != max_distances.Cm2MaxDistance:
            # Write for tolerance of 2cm.
            txtFile.write('Based on 0.02 m tolerance\n')
            txtFile.write('Max Distance: %s m\n' % round(max_distances.Cm2MaxDistance, 2))
            txtFile.write('Records that contain segments that exceed max distance: %s\n' % record_log.Exceeds2CmCount)
            txtFile.write('Percentage of input records that exceed max distance: %s%%\n' % round(((record_log.Exceeds2CmCount / record_log.FeatureCount) * 100 ), 2))
            txtFile.write('\n')

        if max_distances.PrimaryMaxDistance != max_distances.Cm3MaxDistance:
            # Write for tolerance of 3cm.
            txtFile.write('Based on 0.03 m tolerance\n')
            txtFile.write('Max Distance: %s m\n' % round(max_distances.Cm3MaxDistance, 2))
            txtFile.write('Records that contain segments that exceed max distance: %s\n' % record_log.Exceeds3CmCount)
            txtFile.write('Percentage of input records that exceed max distance: %s%%\n' % round(((record_log.Exceeds3CmCount / record_log.FeatureCount) * 100 ), 2))
            txtFile.write('\n')

        if max_distances.PrimaryMaxDistance != max_distances.Cm10MaxDistance:
            # Write for tolerance of 10cm.
            txtFile.write('Based on 0.1 m tolerance\n')
            txtFile.write('Max Distance: %s m\n' % round(max_distances.Cm10MaxDistance, 2))
            txtFile.write('Records that contain segments that exceed max distance: %s\n' % record_log.Exceeds10CmCount)
            txtFile.write('Percentage of input records that exceed max distance: %s%%\n' % round(((record_log.Exceeds10CmCount / record_log.FeatureCount) * 100 ), 2))
            txtFile.write('\n')

        # Write notes.
        txtFile.write('**Note:\n')
        txtFile.write('-All units are converted and calculated in meters.\n')
        txtFile.write('-For a given tolerance value, max distance is the maximum distance\n')
        txtFile.write(' between two vertices above which the route segment is considered as\n')
        txtFile.write(' a segment with sparse vertices.\n')
        txtFile.write('-For the current records, the smallest of the two, XY or Z tolerance\n')
        txtFile.write(' values is considered.\n')
        txtFile.write('-It is recommended that for the densification process, a value smaller\n')
        txtFile.write(' by 10-20% than the max distance should be specified as an acceptable\n')
        txtFile.write(' maximum distance between the vertices.\n')

        txtFile.write('\n')
        txtFile.close()

    # Add message of log file location. 
    arcpy.AddMessage('Log file at %s' % outputFileLog)

##******************************************##
##            Utility Functions             ##
##******************************************##

def GetMaxDistances(in_feature_class, sr):
    
    # Get semi minor axis.
    R = sr.semiMinorAxis;

    # Get Z tolerance in meters.
    Z = sr.ZTolerance
    if sr.VCS is not None:
        unitName = sr.VCS.linearUnitName
        if unitName.lower() != "meter":
            Z = sr.VCS.metersPerUnit

    # Converts a distance in degrees along the equator to the equivalent number of meters.
    XY = sr.XYTolerance * (20015077 / 180 )
 
    # Max Distance = 2 * sqrt(pow(R,2) - pow(R-tolerance,2)).
    maxXY = 2 * (math.sqrt(math.pow(R,2) - math.pow((R-XY),2)))
    maxZ = 2 * (math.sqrt(math.pow(R,2) - math.pow((R-Z),2)))
    maxft = 2 * (math.sqrt(math.pow(R,2) - math.pow((R-0.0003048),2)))
    max1mm = 2 * (math.sqrt(math.pow(R,2) - math.pow((R-0.001),2)))
    max1cm = 2 * (math.sqrt(math.pow(R,2) - math.pow((R-0.01),2)))
    max2cm = 2 * (math.sqrt(math.pow(R,2) - math.pow((R-0.02),2)))
    max3cm = 2 * (math.sqrt(math.pow(R,2) - math.pow((R-0.03),2)))
    max10cm = 2 * (math.sqrt(math.pow(R,2) - math.pow((R-0.1),2)))

    # For the output feature class use the smallest of the XY and Z tolerance. 
    primaryMaxDistance = maxXY if XY < Z else maxZ

    return MaxDistances(primaryMaxDistance, maxXY, maxZ, maxft, max1mm, max1cm, max2cm, max3cm, max10cm)

def ValidateFeatureClass(in_feature_class):

    # Ensure feature exist.
    if arcpy.Exists(in_feature_class):

        # Only features are supported.
        fcPath = arcpy.Describe(in_feature_class).catalogPath
        if 'https:' in fcPath:
            return "Feature layers or layers from a service are not supported"

        # Only feature datasets supported (Feature Class).
        dataType = arcpy.Describe(fcPath).dataType
        if dataType.lower() != "featureclass":
            return "Input is not a valid feature class"

        # Only GCS (unprojected) data supported.
        sp = arcpy.Describe(fcPath).spatialReference
        if sp.type == 'Projected':
            return "Input feature class spatial reference is not in GCS"
    else:
        return "Could not find Feature Class"

def ValidateOutputFeatureClass(path):

    workspaceDir = os.path.dirname(path)
    fileName = path.split('\\')[-1]

    # remove schema prefix if db is sde.
    workspaceType = arcpy.Describe(workspaceDir).workspaceType
    if workspaceType == "RemoteDatabase":
        fileName = fileName.split('.', 1)[-1]

    # make sure out is workspace or valid folder
    dataType = arcpy.Describe(workspaceDir).dataType
    if dataType.lower() != "folder" and dataType.lower() != "workspace":
        return "Output location is not valid"

    # validate if name is valid. If not outfc will be a different name.
    outfc = arcpy.ValidateTableName(fileName, workspaceDir)
    if fileName != outfc:
        return "Output name is not valid"

def CreateOutputFeatureClass(out_path_and_name, in_feature_class, sr):
    
    workspaceDir = os.path.dirname(out_path_and_name)
    fileName = out_path_and_name.split('\\')[-1]

    # Check for m and z awareness.
    desc = arcpy.Describe(in_feature_class)
    hasZ = "ENABLED" if desc.hasZ else "DISABLED"
    hasM = "ENABLED" if desc.hasM else "DISABLED"

    # Create Feature Class.
    arcpy.management.CreateFeatureclass(workspaceDir, fileName, "POLYLINE", in_feature_class, 
                                hasM, hasZ, sr)

    # Return path to the feature.
    return workspaceDir + "\\" + fileName

##******************************************##
##              Named Tuples                ##
##******************************************##

class MaxDistances(NamedTuple):
    PrimaryMaxDistance: float
    XYMaxDistance: float
    ZMaxDistance: float
    FtMaxDistance: float
    Mm1MaxDistance: float
    Cm1MaxDistance: float
    Cm2MaxDistance: float
    Cm3MaxDistance: float
    Cm10MaxDistance: float

class RecordLog(NamedTuple):
    FeatureCount: int
    ExceedsXYCount: int
    ExceedsZCount: int
    ExceedsFtCount: int
    Exceeds1MmCount: int
    Exceeds1CmCount: int
    Exceeds2CmCount: int
    Exceeds3CmCount: int
    Exceeds10CmCount: int

