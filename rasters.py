#
# rasters.py
# Rasters.py contains objects and methods used to open and process raw raster data into a CSV format that can later be
# analyzed by the user using python built-in methods in phenology.py, or using modules written in GNU R.
# 
# Author: Kyle Taylor (kyle.a.taylor@gmail.com) (c) 2012
#

########################################################################################################################
# RELEASED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE, V2.  PLEASE SEE: http://www.gnu.org/licenses/gpl-2.0.txt #
#                                        FOR A FULL TEXT COPY OF THE LICENSE.                                          #
########################################################################################################################

# External Imports

try:
  import operator, sys
  from   osgeo import gdal, gdalnumeric, ogr, osr
  import Image, ImageDraw
  
except ImportError:
  sys.stdout.write("\nWhoopsie Caught from rasters.py\nFailed to import some requisite python libraries.  The OPE python interface requires you have GDAL/OGR and the Python Image Library\ninstalled.  These libraries came pre-packaged with OPE for POSIX and Windows platforms, so if you installed from a standard OPE release,\nyou shouldn't see this message.  See the included install.txt for more information on how to set them up.\n\n")
  sys.exit(-1)

# Local Imports

  from engineCore import dataset
  from phenology  import phenologyCSV_output
  
#
# Raster Image Output
#

class rasterRaw_output:
  def __init__():
    return 0

#
# Raster Input File
#

class rasterRaw_input:
  def __init__(self, path):
    
    ## in-class member variables

    self.validExtensions = [".tif", ".img"]
    self.filepath        = str()
    
    # don't assume the path has been previously checked. Catch any errors and throw an error-quit

    try:
      self.srcArray        = gdalnumeric.LoadFile(path)
      self.srcImage        = gdal.Open(raster)
    except:
      sys.stdout.write("error opening " + path + "... does the file exist / is it a valid image file?\n")
      sys.exit(-1)
    
    self.minIntensity    = float(0)
    self.maxIntensity    = float(0)
    self.meanIntensity   = float(0) 
    
  #
  # pilArrayToGDALnumeric()
  # Converts a Python Imaging Library array to a gdalnumeric image.
  # Re-implemented from code written by J. Lawhead
  #
  
  def pilArrayToGDALnumeric(self, pA):

    gdalArray=gdalnumeric.fromstring(pA.tostring(),'b')
    gdalArray=shape=pA.im.size[1], pA.im.size[0]
    return a
    
  def GDALnumericToPilArray(self, a):
    
    i=Image.fromstring('L',(a.shape[1],a.shape[0]),
            (a.astype('b')).tostring())
    return i

  #
  # geoCoordToPixelLocation()
  # Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
  # the pixel location of a geospatial coordinate 
  # Re-implemented from code written by J. Lawhead
  # 
  
  def geoCoordToPixelLocation(self):
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / xDist)
    return (pixel, line) 

  #
  # buildHistogram()
  # Histogram function for multi-dimensional array.
  # Re-implemented from code written by J. Lawhead
  #
  
  def buildHistogram(self):
    fa = a.flat
    n = gdalnumeric.searchsorted(gdalnumeric.sort(fa), bins)
    n = gdalnumeric.concatenate([n, [len(fa)]])
    hist = n[1:]-n[:-1] 
    return hist
  #
  # stretch()
  # Performs a histogram stretch on a gdalnumeric array image.  
  # Re-implemented from code written by J. Lawhead
  #
  
  def stretch(self):  
    hist = histogram(a)
    im = arrayToImage(a)   
    lut = []
    for b in range(0, len(hist), 256):
      # step size
      step = reduce(operator.add, hist[b:b+256]) / 255
      # create equalization lookup table
      n = 0
      for i in range(256):
        lut.append(n / step)
        n = n + hist[i+b]
    im = im.point(lut)
    return imageToArray(im)
    
  #
  # fileIsRaster()
  # verify that the infile passed to our constructor is actually a raster file.
  #
  
  def fileIsRaster(self, filepath):
    return 0

#
# Mask
#

class mask:
  
  #
  # Constructor __init__
  # 
  
  def __init__(self, maskPath, inputDatasetPath):
    
    ## member variables
    self.inputDataset  = dataset(inputDatasetPath)
    self.outputDataset = list()
    # re-code the dataset class to function as an empty directory that accepts our clipped output rasters.  Right now, 
    # it can only really process input datasets.
    self.clippedOutputDataset = self.inputDataset.inputFileArray[0][:self.inputDataset.inputFileArray[0].rfind("/")] + str("clipped.tif")
    
    if os.path.exists(maskPath):
      self.maskFile = maskPath
    else:
      sys.stdout.write("invalid path specified for mask file: " + maskPath)
      sys.exit(-1)


  #
  # apply()
  # apply a mask file to an input raster file or dataset containing a list of raster files.
  # Re-implemented from code written by J. Lawhead
  #
  
  def apply(self, arg):
    
    # iterate over our input raster dataset, masking each file.
    for raster in self.inputDataset.inputFileArray:

      rasterInput = rasterRaw_input(raster)
      
      # Append a new clip file path for our current raster
      self.outputDataset.append(raster[:raster.rfind(".")]+"clipped"+raster[raster.rfind("."):])

      # Also load as a gdal image to get geotransform 
      # (world file) info

      geoTrans = srcImage.GetGeoTransform()

      # Create an OGR layer from a boundary shapefile
      shapef = ogr.Open("%s.shp" % shp)
      lyr = shapef.GetLayer(shp)
      poly = lyr.GetNextFeature()

      # Convert the layer extent to image pixel coordinates
      minX, maxX, minY, maxY = lyr.GetExtent()
      ulX, ulY = world2Pixel(geoTrans, minX, maxY)
      lrX, lrY = world2Pixel(geoTrans, maxX, minY)

      # Calculate the pixel size of the new image
      pxWidth = int(lrX - ulX)
      pxHeight = int(lrY - ulY)

      clip = srcArray[:, ulY:lrY, ulX:lrX]

      # Create a new geomatrix for the image
      geoTrans = list(geoTrans)
      geoTrans[0] = minX
      geoTrans[3] = maxY

      # Map points to pixels for drawing the 
      # boundary on a blank 8-bit, 
      # black and white, mask image.
      
      points = []
      pixels = []
      geom = poly.GetGeometryRef()
      pts = geom.GetGeometryRef(0)
      for p in range(pts.GetPointCount()):
        points.append((pts.GetX(p), pts.GetY(p)))
      for p in points:
        pixels.append(world2Pixel(geoTrans, p[0], p[1]))
      rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
      rasterize = ImageDraw.Draw(rasterPoly)
      rasterize.polygon(pixels, 0)
      mask = imageToArray(rasterPoly)   

      # Clip the image using the mask
      clip = gdalnumeric.choose(mask, \
          (clip, 0)).astype(gdalnumeric.uint8)

      # This image has 3 bands so we stretch each one to make them
      # visually brighter
      for i in range(3):
        clip[i,:,:] = stretch(clip[i,:,:])

      # Save ndvi as tiff
      gdalnumeric.SaveArray(clip, "%s.tif" % output, \
          format="GTiff", prototype=raster)

      # Save ndvi as an 8-bit jpeg for an easy, quick preview
      clip = clip.astype(gdalnumeric.uint8)
      gdalnumeric.SaveArray(clip, "%s.jpg" % output, format="JPEG")
      
#
# Core Methods for Analysis of Raster Data
#

class rasterCalculator:
  def blah():
    return 0

