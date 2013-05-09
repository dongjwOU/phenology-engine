#
# engineCore.py
#
# Core implementation for the Open Phenology Engine.  Contains the default console interface to engine components
# for downloading, masking, and processing datasets.
#
# Author: Kyle Taylor (kyle.a.taylor@gmail.com) (c) 2012
#

########################################################################################################################
# RELEASED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE, V2.  PLEASE SEE: http://www.gnu.org/licenses/gpl-2.0.txt #
#                                        FOR A FULL TEXT COPY OF THE LICENSE.                                          #
########################################################################################################################

import sys,os,subprocess
import math

#
# local imports
#

from phenology import *
from rasters   import *
from net       import *

#
# errorReporting()
# methods for handling console feedback, verbosity, and eventually crashes
#

class errorReporter:
  def __init__(self, *args):

    try:
      self.verbosity = int(args[0])
    except IndexError:
      self.verbosity = int(0)

  def setVerbosity(self, arg):
    self.verbosity = int(arg)
  def write(self, message):
    if(self.verbosity > 0):
      sys.stdout.write(message)

#
# csvDataHandler()
# csvDataHandler contains methods used to manipulate a CSV input/output file.  Think of it as a simpler,
# straight-foward interface to phenologyCSV_input / phenologyCSV_output classes imported from phenology
#

class csvDataHandler:

  # in-class methods
  def __init__(self, args):

    self.mode = str()

    # let's attempt to interact with whatever the user passed as if it were
    # a valid _csv object.  If it isn't a valid object, we can throw an exception

    if(os.path.exists(args.path) == False):
      print " couldn't find path: " + str(args)
      sys.exit(-1)

    else:

      # perhaps the user passed a CSV object?
      # let's try and work with it to determine
      # whether its a phenology input/output csv

      try:
        # this is a phenologyCSV_infile object...
        if(args.infile):
          try:
            self.raw = args
            self.setReadMode(self.raw)
            self.mode = "r"
          except Exception, e:
            print e
      except AttributeError:
        # this is a phenologyCSV_outfile object...
        if(args.outfile):
          try:
            self.raw = args
            self.setWriteMode(self.raw)
            self.mode = "w"
          except Exception, e:
            print e
      except Exception, e:
        print e
        sys.exit(-1)

  #
  # setWriteMode((self, args)
  # Assumes an output file has already been created / opened by phenologyCSV_input.  This just hijacks an instance
  # of that file and makes it available in setWriteMode( mode through csvDataHandler as csvDataHandler.File
  #

  def setWriteMode(self, args):
    self.raw = args
    self.File = args.outfile

  #
  # setReadMode(self, args)
  # Opens a phenology csv file in read mode and verifies the header structure to ensure that we are working with a properly formated CSV file.
  # a file object is made avaiable as csvDataHandler.File
  #

  def setReadMode(self, args):

    if(len(str(args)) > 0):

      # If we haven't defined a file yet, let's do it now.
      # The default action is for the first argument passed to setReadMode() to be
      # an initialized _csv file object.

      if(len(self.raw.path) > 0):

        if(args.infile):
          try:
            self.raw  = args
            self.File = args.infile
          except Exception, e:
            print e
        else:
          print "error.  no valid input file specified to read."
          sys.exit(-1)

    # verify the header signature to make sure that this is a CSV file
    # we know how to work with.

    line = self.File.readline()

    if(line.strip() != self.raw.defaultHeader.strip()):
        print "invalid header found in input CSV file"
        sys.exit(-1)

    # Now process the core data of the file

    while line:

      elements = line.split(",")

      if elements[0] != "Start_Day": # Make sure we aren't processing header data

        self.raw.start_days.append(int(elements[0]))
        self.raw.mean_intensities.append(float(elements[3]))

        # correct for ranges that overlap into the next year...

        if ( (int(elements[0]) >= 358) & (int(elements[1]) < 10) ):
          self.raw.end_days.append(int(elements[0]) + 7)
          self.raw.median_days.append(int(elements[0]) + int(elements[1])  / 2 )
        else:
          self.raw.end_days.append(int(elements[1]))
          self.raw.median_days.append(int(elements[2]))

      line = self.File.readline()

#
# Dataset class contains variables and methods associated with processing input data
# for the phenology engine.
#

class dataset:

  # inputFileArray is instantiated as an empty list, but it will eventually hold a list
  # of csvDataHandler::phenologyCSV_input objects or a list of rasterDataHandler::raster_input
  # objects.

  working_dir = "N"
  script_path = "N"

  #
  # __init__() Constructor
  #
  # accepts a character string path to either a directory or single file.  If it is a directory path, it will
  # crawl the contents of the directory and append individual files to a list (inputFileArray).  Do not assume
  # that the input path string is valid.  Error-quit on invalid path.
  
  def __init__(self, args):

    self.inputFileArray = list()

    # did the user pass a single file?
    if(os.path.isfile(args)):
      inputCSV = phenologyCSV_input(str(args))
      #self.inputFileArray.append(inputCSV)

    # did the user pass a whole directory path?
    elif(os.path.exists(args)):
      for root, dirs, files in os.walk(args):
        if self.working_dir == "N": self.working_dir = root
        for f in files:
          inputCSV = phenologyCSV_input(root+"/"+f)
          self.inputFileArray.append(inputCSV)

    else:
      print " path not found: ", args
      sys.exit(-1)

  #
  # parse()
  # extract whatever data was stored in inputFileArray into our application space
  #

  def parse(self):
    # iterate over the elements of the file array.
    for f in self.inputFileArray:
      # is this a CSV input file?
      if f.path[len(f.path)-4:] == ".csv":
        f.extract()
      # is this a raster input file?

#
# engineParameters
# Contains initial run-time settings that determine how the engine operates
#

class engineParameters:
  def __init__(self):

    ## in-class member variables

    self.OPERATION_MODE           = "REPORT"  # Let's report by default, as I can't get GDAL to mask yet...
    self.FOUND_RASTER_INPUT_DATA  = str()
    self.FOUND_CSV_INPUT_DATA     = str()
    self.USING_SNOW_DATA          = str()
    self.USING_CUSTOM_FUZZ        = "N"  # Are we using a custom fuzz factor?
    self.USING_OUT_CSV_OPTION     = "N"  # Is engine output directed to a CSV file?
    self.USING_INPUT_DATASET      = "N"
    self.USING_SNOW_DATA          = "N"
    self.USING_UNIMODAL_EXTRACT   = "N"
    self.FUZZ_MAX                 = 1500  # Set the default FUZZ_FACTOR value
    self.FUZZ_FACTOR              = 0.1
    self.VERBOSITY                = 0
    self.VERSION                  = 0.1

    self.inputDataset             = str()  # temporarily initiated as an empty string.  Will eventually become a dataset()
    self.snowDataset              = str()  #

  #
  # validate()
  # Runs through a series of checks to ensure that all appropriate options are enables and runtime variables are set
  # before initiating the phenology engine.  Returns true/false

  def validate(self):

    # chiefly, do we have an input dataset to work with?

    if self.USING_INPUT_DATASET == "N" and self.USING_SNOW_DATA == "N":
      sys.stdout.write("no input dataset defined...\n")
      sys.exit(-1)

    # fetch our operation-mode and determine if it is consistent with the input data the user provided.

    if self.OPERATION_MODE == "REPORT":
     if self.USING_INPUT_DATASET == "Y":
      for f in self.inputDataset.inputFileArray:
        if f.path[len(str(f.path))-4:] != ".csv":
          print "\n fatal error: input data " + f + " is not valid for processing a phenology report but is in the directory structure you passed at runtime.\n"
          sys.exit(-1)
     elif self.USING_SNOW_DATA == "Y":
      for f in self.snowDataset.inputFileArray:
        if f.path[len(str(f.path))-4:] != ".csv":
          print "\n fatal error: input data " + f + " is not valid for processing a phenology report but is in the directory structure you passed at runtime.\n"
          sys.exit(-1)


  def printUsage(self):
    print "\nPHENOLOGY ENGINE" + " (v" + str(self.VERSION) + " alpha) - http://phenology-engine.googlecode.com"
    print "Kyle Taylor (c) 2012 - Released Under the Terms of the GNU General Public License (v2)"
    print "-----------------------------------------------------------------------------------------------------------------------------------------"
    print "OPTIONS: --dataset <path>     : specify an input file or directory path containing input files."
    print "         --operation <string> : BUILD, EXTRACT, or REPORT.  Without this option, engine will default to BUILD.  See README for more info."
    print "         --module <path>      : accept an R script module through which we pass engine output.  See README for usage."
    print "         --fuzz <number>      : set the fuzz factor used for finding days corresponding to raster intensities."
    print "         --fetch <string>     : attempt to fetch raw satellite data from USGS FTP servers to BUILD with.  See README for query syntax."
    print "         --snowdata <path>    : use supplemental snow data at <path> to improve the predictions of our days."
    print "         --mask <filename>    : mask input data with <filename>.  This option is required for BUILD operations."
    print "         --output <filename>  : dump information generated to <filename>."
    print "         --verbose            : print a lot of output to the console to keep track of what the engine is doing."
    print "         --unimodal-repair    : attempt to extract a single unimodal distribution from source data that may have many noisy peaks."
    print "\nEXAMPLE:  " + sys.argv[0] + " --dataset ndvi_raster_data/ --operation BUILD --mask plot_raster_mask.shp --output report.csv\n"
    sys.exit(-1)


#
# __ main() __
#

def main():

 # Create an instance of our engine runtime settings
 RuntimeSettings = engineParameters()
 console         = errorReporter()

 if len(sys.argv) < 2:
   RuntimeSettings.printUsage()

 #
 # Process runtime arguments and options passed by the user
 #

 for i in xrange(len(sys.argv)):

   if sys.argv[i] == "--help":
     RuntimeSettings.printUsage()

   elif sys.argv[i] == "--verbose":
     RuntimeSettings.VERBOSITY = 1
     console.setVerbosity(1)

   elif sys.argv[i] == "--unimodal-repair":
     RuntimeSettings.USING_UNIMODAL_EXTRACT = "Y"

   elif sys.argv[i] == "--fuzz":
     RuntimeSettings.MAX_FUZZ = float(sys.argv[i+1])
     RuntimeSettings.USING_CUSTOM_FUZZ = "Y"

   elif sys.argv[i] == "--output":
     # sanity check the next argument
     if(i+1 >= len(sys.argv)):
       print "--output requires a filename."
       sys.exit(-1)

     # assign runtime settings and output file
     RuntimeSettings.USING_OUT_CSV_OPTION = "Y"
     outCSV = phenologyCSV_output(str(sys.argv[i+1]))
     output = csvDataHandler(outCSV)

   # --dataset will accept a PATH string to either a file or a directory
   # if the user passes a directory path, we are going to need to scan it for its contents

   elif sys.argv[i] == "--dataset":
     # sanity check the next argument
     if(i+1 >= len(sys.argv)):
       print "--dataset requires a path"
       sys.exit(-1)

     # assign runtime settings and input data
     RuntimeSettings.USING_INPUT_DATASET = "Y"
     RuntimeSettings.inputDataset = dataset(sys.argv[i+1])

   elif sys.argv[i] == "--snowdata":
     RuntimeSettings.USING_SNOW_DATA = "Y"
     RuntimeSettings.snowDataset = dataset(sys.argv[i+1])

 ## Validate Runtime Settings

 RuntimeSettings.validate()

 ## Fetch our operational mode and act accordingly

 if RuntimeSettings.OPERATION_MODE == "BUILD":   # parse input CSV data into our application space
   RuntimeSettings.inputDataset.parse()
   print "BUILD not yet implemented... we're getting there."
   sys.exit(-1)
   # when finished, set our operation mode to the next step in the process.
   RuntimeSettings.OPERATION_MODE = "EXTRACT"
 if RuntimeSettings.OPERATION_MODE == "EXTRACT":
   print "EXTRACT not yet implemented... we're getting there."
   sys.exit(-1)
   RuntimeSettings.OPERATION_MODE = "REPORT"
 if RuntimeSettings.OPERATION_MODE == "REPORT":

   # process snow data
   if RuntimeSettings.USING_SNOW_DATA == "Y":

    for infile in RuntimeSettings.snowDataset.inputFileArray:
      console.write("(snow="+str(infile.path[infile.path.rfind("/")+1:])+")")

      snowSeason = phenologyCalculator(infile)

      # normalize the snow distribution so our algorithms can work with the dataset
      snowSeason.correctDistFor(snowSeason.minVal)

      # find the start/end days of snow season
      snowSeasonInitiation  = snowSeason.retFirstHalfSeason("snow")
      snowSeasonTermination = snowSeason.retSecondHalfSeason("snow")

      snowSeasonInitiation.reverse()

      # we use big co-factors here (100), because the jump between snow cover / uncover is usually a huge jump.
      snowInitDay = snowSeason.pettorelliRunningAvgMethod(snowSeasonInitiation, 100)
      snowTermDay = snowSeason.pettorelliRunningAvgMethod(snowSeasonTermination, 100)

      console.write("(snow free day: "+str(snowInitDay)+")"+"(snow cover day: "+str(snowTermDay)+");")

      if RuntimeSettings.USING_INPUT_DATASET == "N" and RuntimeSettings.USING_OUT_CSV_OPTION == "Y":
        output.raw.outfile.write(str(infile.year)+","+str(infile.unit)+","+str(infile.plot)+",,,,,,,,,,,"+str(snowInitDay)+","+str(snowTermDay)+",,,,\n")

  # process core AVHRR/MODIS vegetation dataset

   if RuntimeSettings.USING_INPUT_DATASET == "Y":
    for infile in RuntimeSettings.inputDataset.inputFileArray:

      console.write("(in="+str(infile.path[infile.path.rfind("/")+1:])+")")

      season = phenologyCalculator(infile)

      # determine how noisey our data are?

      # divide our distribution into initiation and termination periods,
      # (mostly) symmetric about seasonal max day.
      seasonInitiation  = season.retFirstHalfSeason("veg")
      seasonTermination = season.retSecondHalfSeason("veg")

      # Find green-up event
      runAvgGreenUpDay = season.pettorelliRunningAvgMethod(seasonInitiation, 1.5)
      staticGreenUpDay = season.findVal(seasonInitiation, float(season.maxVal)*0.2, RuntimeSettings.FUZZ_FACTOR, RuntimeSettings.FUZZ_MAX)

      console.write("(run. avg. green-up day: " + str(runAvgGreenUpDay) + ")")
      console.write("(p=0.2 green-up day: " + str(staticGreenUpDay) + ")")

      # Find senescence event
      season.seasonTermination.reverse() # reverse the late season curve so it matches the logic of the early season curve
      runAvgSenescenceDay = season.pettorelliRunningAvgMethod(seasonTermination, 1.5)
      staticSenescenceDay = season.findVal(seasonTermination, float(season.maxVal*0.3), RuntimeSettings.FUZZ_FACTOR, RuntimeSettings.FUZZ_MAX)

      console.write("(run. avg. senescence day: " + str(runAvgSenescenceDay) + ")")
      console.write("(p=0.3 senescence day: " + str(staticSenescenceDay) + ")")

      # perform a unimodal extraction of the source data if the user requested it.  This will hopefully prevent our halfmax day values
      # from getting cluttered when analyzing noisy distributions

      if RuntimeSettings.USING_UNIMODAL_EXTRACT == "Y":
        # create a data filter and then re-process our season initiation and termination distributions
        distFilter = phenologyDataFilter()
        seasonInitiation  = distFilter.extractUnimodal(seasonInitiation, runAvgGreenUpDay, season.maxDay)
        seasonTermination = distFilter.extractUnimodal(seasonTermination, season.maxDay, runAvgSenescenceDay)

      # Find half max days
      halfMaxGreenUp = season.findVal(seasonInitiation, float(season.maxVal*0.5), RuntimeSettings.FUZZ_FACTOR, RuntimeSettings.FUZZ_MAX)
      halfMaxSenescence = season.findVal(seasonTermination, float(season.maxVal*0.5), RuntimeSettings.FUZZ_FACTOR, RuntimeSettings.FUZZ_MAX)

      console.write("(halfmax green-up day: " + str(halfMaxGreenUp)+")")
      console.write("(halfmax senescence day: " + str(halfMaxSenescence)+")")

      # Calculate Half-Max LOS
      halfMaxLOS = halfMaxSenescence - halfMaxGreenUp

      # Do some sanity checks on our results

      console.write("(halfmax LOS: " + str(halfMaxLOS)+")")

      console.write(";")

      # Write to an Output CSV report
      # "Year, Unit, Plot, Peak_Intensity, p50_Intensity, NDVI_Mean_Intensity, Run_Avg_Green_up_Day, P=0.2_Green_up_Day, Mid_Slope_of_Increase_Day, Peak_Intensity_Day, Mid_Slope_of_Decrease_Day, Run_Avg_Senescence_Day, P=0.3_Senescence_Day, Last_P90_Snowcover_Start, Last_P90_Snowcover_End, Green_Up_LOS, Half_Max_LOS, Fuzz_Factor, Sanity_Check"

      if RuntimeSettings.USING_OUT_CSV_OPTION == "Y":
        # with snow data
        if RuntimeSettings.USING_INPUT_DATASET == "Y" and RuntimeSettings.USING_SNOW_DATA == "Y":
          output.raw.outfile.write(str(infile.year)+","+str(infile.unit)+","+str(infile.plot)+","+str(season.maxVal)+","+str(float(season.maxVal)/2)+","+str(season.findMean())+","+str(runAvgGreenUpDay)+","+str(staticGreenUpDay)+","+str(halfMaxGreenUp)+","+str(season.maxDay)+","+str(halfMaxSenescence)+","+str(runAvgSenescenceDay)+","+str(staticSenescenceDay)+","+str(snowInitDay)+","+str(snowTermDay)+","+str(int(runAvgSenescenceDay-runAvgGreenUpDay))+","+str(halfMaxLOS)+","+str(RuntimeSettings.FUZZ_FACTOR)+"\n")
        # without snow data
        if RuntimeSettings.USING_INPUT_DATASET == "Y" and RuntimeSettings.USING_SNOW_DATA == "N":
          output.raw.outfile.write(str(infile.year)+","+str(infile.unit)+","+str(infile.plot)+","+str(season.maxVal)+","+str(float(season.maxVal)/2)+","+str(season.findMean())+","+str(runAvgGreenUpDay)+","+str(staticGreenUpDay)+","+str(halfMaxGreenUp)+","+str(season.maxDay)+","+str(halfMaxSenescence)+","+str(runAvgSenescenceDay)+","+str(staticSenescenceDay)+","+str(" ")+","+str(" ")+","+str(int(runAvgSenescenceDay-runAvgGreenUpDay))+","+str(halfMaxLOS)+","+str(RuntimeSettings.FUZZ_FACTOR)+"\n")

if __name__ == "__main__":
    main()
    sys.stdout.write("\n(quit)\n")
