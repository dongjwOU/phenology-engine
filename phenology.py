#
# phenology.py
# Objects and methods associated with processing comma-separated file data previously compiled from raster data (see: rasters.py).
# Phenology.py's algorithms are the core of our analysis engine, as they are what are employed to analyse our source datasets.
# 
# Author: Kyle Taylor (kyle.a.taylor@gmail.com) (c) 2012
#

########################################################################################################################
# RELEASED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE, V2.  PLEASE SEE: http://www.gnu.org/licenses/gpl-2.0.txt #
#                                        FOR A FULL TEXT COPY OF THE LICENSE.                                          #
########################################################################################################################


import sys, os
import math

#
# distribution()
# an object to define intensities and their corresponding day-of-year without having to return
# cumbersome tuples from functions.
#

class distribution:
  def __init__(self):
    self.days        = list()
    self.intensities = list()
  def reverse(self):
    self.days.reverse()
    self.intensities.reverse()

#
# Phenology CSV Input Data
#

class phenologyCSV_input:
  def __init__(self, args):

    self.path               = str()
    self.defaultHeader      = "Start_Day, End_Day, Median_Day, Mean, Stdev"

    self.start_days         = list()
    self.end_days           = list()
    self.median_days        = list()
    self.raster_intensities = list()

    self.year               = int() # gleaned from filename
    self.plot               = str() #
    self.unit               = str() #

    # First argument should be a string containing the path to a source csv
    # Verify that the path exists.

    if( len(args) > 0 ):
      if(os.path.exists(args)):
        self.path = str(args)
      else:
        print " couldn't find path: " + str(args)
        sys.exit(-1)
    else:
      print "no args passed at object initiation for PhenologyCSV_input::__init__()"
      sys.exit(-1)
    try:
      self.infile = open(self.path, "r")
      self.verifyHeader()

    except IOError, e:
      print e.message
      sys.exit(-1)

    # Try and establish the year and sampleLabel from the infile name

    if(self.path.rfind("/") != -1): # directories will have "/" in their path
      fileNameLabel = self.path[self.path.rfind("/")+1:]
    else: # single files in the CWD will have no "/" in their path
      fileNameLabel = self.path

    # find year
    if fileNameLabel.find("-") == 4: # The first 4 characters of an infile should always be the year, followed by "-"
      self.year = int(fileNameLabel[:4])
    else:
      sys.stdout.write("(couldn't determine year from infile: " + fileNameLabel + ")")
      sys.exit(-1)

    # find sample label

    if(self.path.find("unit") != -1): # are we dealing with sample units?
      self.unit = fileNameLabel[fileNameLabel.find("unit")+5:]
      self.unit = self.unit[:self.unit.find("-")]

    if(self.path.find("plot") != -1): # are we dealing with sample plots?
      self.plot = fileNameLabel[fileNameLabel.find("plot")+5:]
      self.plot = self.plot[:self.plot.find("-")]

    # Read-in data from infile into variable space of our application

    for ln in self.infile:

      elements = ln.split(",")

      if elements[0] != "Start_Day": # Make sure we aren't processing header data

        self.start_days.append(int(elements[0]))
        self.raster_intensities.append(float(elements[3]))

        # correct for ranges that overlap into the next year...

        if ( (int(elements[0]) >= 358) & (int(elements[1]) < 10) ):
          self.end_days.append(int(elements[0]) + 7)
          self.median_days.append(int(elements[0]) + int(elements[1])  / 2 )
        else:
          self.end_days.append(int(elements[1]))
          self.median_days.append(int(elements[2]))

  #
  # verifyHeader()
  # Read the first line of an open infile, strip its white space, and match it against our defaultHeader
  #

  def verifyHeader(self):
   self.infile.seek(0)

   if(self.infile.readline().strip() != self.defaultHeader):
     return False
   else:
     return True
#
# Phenology CSV Output Data
#

class phenologyCSV_output:
  def __init__(self, args):
    ## in-class member variables

    self.outfile       = str()
    self.path          = str()
    self.defaultHeader = "Year, Unit, Plot, Peak_Intensity, p50_Intensity, NDVI_Mean_Intensity, Run_Avg_Green_up_Day, P=0.2_Green_up_Day, Mid_Slope_of_Increase_Day, Peak_Intensity_Day, Mid_Slope_of_Decrease_Day, Run_Avg_Senescence_Day, P=0.3_Senescence_Day, Last_P90_Snowcover_Start, Last_P90_Snowcover_End, Run_Avg_LOS, Half_Max_LOS, Fuzz_Factor, Sanity_Check"

    # First argument should be a string containing the path to a source csv

    if( len(args) > 0 ):

      # append to existing output file?
      if(os.path.exists(args)):

        self.path = str(args)
        self.outfile = open(self.path, "a")
        sys.stdout.write("(out=" + self.path + ", a);")

      # create a new file and write a default header to it.
      else:

        self.path = str(args)

        try:
          self.outfile = open(self.path, "w")
          self.outfile.write(self.defaultHeader + "\n")
          sys.stdout.write("(out=" + self.path + ",w)(header);")
        except Exception, e:
          print e.message
          sys.exit(-1)
    else:
      print "no args passed at object initiation for phenologyCSV_output::__init__()"
      sys.exit(-1)

#
# Filters
# core methods for filtering distributions from processed (but perhaps noisey) raster input data
#

class phenologyDataFilter:
  def __init__(self):
    return None
  def findNoiseLevel(self):
    return 0
  def normalityTest(self):
    return False

  #
  #
  # extract a single unimodal distribution from a noisy dataset that may have multiple peaks, using previously processed
  # green-up/senescence day data.
  #

  def extractUnimodal(self, seasonSample, greenUpDay, senescenceDay):

    returnThis = distribution()

    # iterate over a distribution and extract all data that fall between green-up and senescence day
    for i in xrange(len(seasonSample.days)):
      if seasonSample.days[i] >= greenUpDay and seasonSample.days[i] <= senescenceDay:
        returnThis.intensities.append(seasonSample.intensities[i])
        returnThis.days.append(seasonSample.days[i])

    return returnThis

  def gausianSmooth(self):
    print "blah"
  def savitzkyGolaySmooth(self):
    print "blah"
  def logisticSmooth(self):
    print "blah"

#
# Phenology Calculator()
# core methods for analysis of processed raster distributions
#

class phenologyCalculator:
  def __init__(self, arg):

    try:

      # the user can pass a phenologyCalculator() object if they choose.  Attempt to work with such an object.
      # If they passed a file() object, we can deal with it by catching an AttributeError and acting accordingly

      if len(arg.sampleDist) > 1:
        self = arg

    except AttributeError:

      self.sampleDist = list() # raster values for whole season
      self.sampleDays = list() # corresponding days for whole season

      self.seasonInitiation = distribution()
      self.seasonTermination = distribution()

      self.minVal   = float(0)
      self.minDay   = int(0)
      self.maxVal   = float(0)
      self.maxDay   = int(0)

      # make sure that we are getting good input data
      if len(arg.raster_intensities) < 1 or len(arg.median_days) < 1:
        print "(no valid in data found from input file: " + str(arg.path) + ")"
      else:
        try:
          for i in xrange(len(arg.raster_intensities)):
            self.sampleDist.append(float(arg.raster_intensities[i]))
            self.sampleDays.append(int(arg.median_days[i]))
        except AttributeError:
          print "(phenologyCalculator() fatal : " + str(arg.path) + "  constructor didn't contain sample data)"
          sys.exit(-1)

    # improve the curve of our distribution by subtracting its minimum value
    self.correctDistFor(min(self.sampleDist))

    # calculate some basic information about our distribution that our member methods depend on
    self.findMin()
    self.findMax()

  #
  # findMin()
  # Iterate over a user-provided sampleDist and look for the maximum/minimum value and corresponding day of that value
  #

  def findMin(self):

    # sanity check
    if len(self.sampleDist) < 1 or len(self.sampleDays) < 1:
      return [-1, -1]

    # set the first day of the distribution as the default minimum value
    self.minVal = float(self.sampleDist[0])
    self.minDay = int(self.sampleDays[0])

    for i in xrange(len(self.sampleDist)):
      if float(self.sampleDist[i]) < float(self.minVal):
        self.minVal = float(self.sampleDist[i])
        self.minDay = int(self.sampleDays[i])
        returnThis = [self.minVal, self.minDay]

    return [self.minVal,self.minDay]

  #
  # findMax()
  #
  #

  def findMax(self):

    # sanity check
    if len(self.sampleDist) < 1 and len(self.sampleDays) < 1:
      return [-1, -1]

    # set the first day of the distribution as the default max value
    self.maxVal = float(self.sampleDist[0])
    self.maxDay = int(self.sampleDays[0])

    for i in xrange(len(self.sampleDist)):

      if float(self.sampleDist[i]) > float(self.maxVal):
        self.maxVal = float(self.sampleDist[i])
        self.maxDay = int(self.sampleDays[i])
        returnThis = [self.maxVal, self.maxDay]

    return [self.maxVal, self.maxDay]

  #
  # findVal()
  # return the first index of the sampleDist containing "value"
  # (left-sided search function)
  #

  def findVal(self, sampleDist, value, FUZZ_FACTOR, FUZZ_MAX):

    returnVal = -1  # by default, return an invalid sampleDist position on error so we can catch it
    run = True

    while(run):
      if len(sampleDist.intensities) > 1:
        if FUZZ_FACTOR < FUZZ_MAX:

          for i in xrange(len(sampleDist.intensities)):
            #sys.stdout.write("(i: " + str(i) + ")")
            #sys.stdout.write("(search:"+str(sampleDist.intensities[i])+" | fuzz: " + str(FUZZ_FACTOR) + ");")
            #sys.stdout.write("(cond: " + str(math.fabs(float(sampleDist.intensities[i]) - float(value))) + ");")

            if float(sampleDist.intensities[i]) == float(value): # exact match?
              return sampleDist.days[i]
            elif ( math.fabs(float(sampleDist.intensities[i]) - float(value)) < FUZZ_FACTOR ): # fuzzy match
              return sampleDist.days[i]

          # clearly we didn't find it.  Let's play with the FUZZ_FACTOR and try again with recursion...
          FUZZ_FACTOR = float(FUZZ_FACTOR) + float(0.1)

        else:
          sys.stdout.write("(max fuzz reached)")
          run=False

      else:
        sys.stdout.write("(sample size <= 1)")
        run=False

    return returnVal

  #
  # return the last index of the sampleDist containing "value"
  # (right-sided search function)
  #

  def rfindVal(self, sampleDist, value):

    global FUZZ_FACTOR

    lastPosFound = -1

    for i in xrange(len(sampleDist)):
        if float(sampleDist[i]) == float(value):
          lastPosFound = i
        elif ( math.fabs(float(sampleDist[i]) - float(value)) < FUZZ_FACTOR ): # fuzzy match... within 0.9?
          lastPosFound = i

    if lastPosFound == -1 and ( FUZZ_FACTOR < 20 or USING_CUSTOM_FUZZ == "Y"): # A fuzz-factor of 15 should never be reached, but we'll use it as a cut-off.
      # clearly we didn't find it.  Let's play with the FUZZ_FACTOR and try again with recursion...
      FUZZ_FACTOR = FUZZ_FACTOR + 0.05
      lastPosFound = rfindValPos(sampleDist, value)

    elif lastPosFound > -1 and FUZZ_FACTOR < 5:
      # we found an occurance, but perhaps the FUZZ_FACTOR was to narrow and we are missing a value close to the end of the sample dist.
      # careful not to over-fuzz by setting this max FUZZ_FACTOR well above 1.
      FUZZ_FACTOR = FUZZ_FACTOR + 0.01
      try:
        lastPosFound = rfindValPos(sampleDist, value)
      except RuntimeError:
        print "--[ERROR]---------------"
        try:
          print "For Year: " + info[0] + ", Unit: " + info[2]
        except:
          print " ..."
        print "... max recursion depth reached without finding a value in rfindValuePos."

    return lastPosFound

  #
  # findMean()
  # Iterate over a sample distribution and calculate the mean.  We could just use numpy, but its really kind of a waste
  # adding an import for something so simple.
  #

  def findMean(self):
    return sum(self.sampleDist) / len(self.sampleDist)

  #
  # retFirstHalfSeason()
  # Iterate over the season and build a list of all values less than the maxvalue

  def retFirstHalfSeason(self, sampleType):

    # most veg data is normally distributed.  MODIS snow data, however,
    # will be a normal distribution flipped on its head.  Account for this.

    for i in xrange(len(self.sampleDays)):
      if sampleType == "veg":
        # simply divide the sample distribution in half and return the first half.  This is a simple approach, but won't work
        # for ecosystems that experience more than one greening season a year, as in many areas of the tropics.  A more elegant
        # system is needed.
        if i <= (len(self.sampleDays)/2):
          self.seasonInitiation.intensities.append(float(self.sampleDist[i]))
          self.seasonInitiation.days.append(int(self.sampleDays[i]))
        # assume that the max value is the middle of a season.  Return everything before the position of that max value.
        #if self.sampleDist[i] <= self.maxVal:  # if we haven't reached the last position of max intensity in the array
          #self.seasonInitiation.intensities.append(float(self.sampleDist[i]))
          #self.seasonInitiation.days.append(int(self.sampleDays[i]))
        else:
          break
      elif sampleType == "snow":
        # simply divide the distribution in half and return the first half.
        if i <= (len(self.sampleDays)/2):
          self.seasonInitiation.intensities.append(float(self.sampleDist[i]))
          self.seasonInitiation.days.append(int(self.sampleDays[i]))
        # assume that the min value is the middle of a season.  Return everything before the position of that min value.
        #if self.sampleDist[i] <= self.minVal:  # if we haven't reached the last position of max intensity in the array
          #self.seasonInitiation.intensities.append(float(self.sampleDist[i]))
          #self.seasonInitiation.days.append(int(self.sampleDays[i]))
        else:
          break

    return self.seasonInitiation

  #
  # retSecondHalfSeason()
  #

  def retSecondHalfSeason(self, sampleType):

    # most veg data is normally distributed.  MODIS snow data, however,
    # will be a normal distribution flipped on its head.  Account for this.

    for i in xrange(len(self.sampleDays)):

      if sampleType == "veg":
        if i >= (len(self.sampleDays)/2):
          self.seasonTermination.intensities.append(float(self.sampleDist[i]))
          self.seasonTermination.days.append(int(self.sampleDays[i]))
        #if self.sampleDays[i] >= self.maxDay:  # if we haven't reached the last position of max intensity in the array
          #self.seasonTermination.intensities.append(float(self.sampleDist[i]))
          #self.seasonTermination.days.append(int(self.sampleDays[i]))
      elif sampleType == "snow":
        # simply cut the distribution in half
        if i >= (len(self.sampleDays)/2):
          self.seasonTermination.intensities.append(float(self.sampleDist[i]))
          self.seasonTermination.days.append(int(self.sampleDays[i]))
        # seek out the min day as the distribution divider
        #if self.sampleDays[i] >= self.minDay:  # if we haven't reached the last position of max intensity in the array
          #self.seasonTermination.intensities.append(float(self.sampleDist[i]))
          #self.seasonTermination.days.append(int(self.sampleDays[i]))

    return self.seasonTermination

  #
  # correct every value of a sample distribution for value, then return a copy of the sampleDist
  #

  def correctDistFor(self, value):
    for i in xrange(len(self.sampleDist)):
      self.sampleDist[i] = self.sampleDist[i] - value

    return self.sampleDist

  #
  # convertSnowDistToNormalDist()
  # if the working distribution is MODIS snow data, it will not be normally distributed.  It will be an inverted normal distribution
  # with the minimum intensity somewhere around day 185 and two peaks at the start and end of the season.  This is experimental.  There is
  # an existing mechanism in retFirstHalfSeason / retSecondHalfSeason that will account for snow data without fundamentally changing the raw
  # NDVI values.

  def convertSnowDistToNormalDist(self):

    # use a copy of the original distribution here.  We are going to significantly change the original recorded values
    # for raster intensities here, and we don't want to potentially mess-up upstream methods that depend on the true raster
    # intensity values in their analysis

    normalizedSnowDist = self.sampleDist

    for i in xrange(len(self.sampleDist)):
      if self.sampleDist[i] > 0:
        normalizedSnowDist[i] = (1 / self.sampleDist[i]) * self.maxVal
      # for all values that equal zero (occurs when a distribution is corrected for its mimimum value
      else:
        normalizedSnowDist[i] = self.maxVal

    return normalizedSnowDist

  #
  # pettorelli_runningAvg()
  # The Pettorelli method involves taking a running average from the start of a distribution and
  # accepting the first day of year where raster intensity jumps +2SD from the running average as
  # the first day of veg change.  This could be a green-up event, or it could be a senescence day if
  # you are coming from the opposite side of the distribution.  It is based loosely on notes in the
  # literature by Pettorelli et al. We implement it here to find green-up and senescence days in distributions
  #

  def pettorelliRunningAvgMethod(self, halfSeason, cofactor):

    while(True):
      runningAverage = 0

      for i in xrange(len(halfSeason.intensities)):
        # zeros are only significant in our distribution at its start.  Noisy zeros in the middle of a
        # valid raster dataset can really throw-off our calculations.  The condition below should fix this.
        if halfSeason.intensities[i] != 0 or (i/len(halfSeason.intensities)) < 0.3 :
          # reset our SD
          standardDeviation = 0

          # calculate the running average
          runningAverage = (runningAverage + halfSeason.intensities[i]) / (i+1)

          # calculate the new standard deviation
          for j in xrange(i): # sum the list of elements - mean
            standardDeviation = standardDeviation + (halfSeason.intensities[j]-runningAverage)
            if i > 1:
              standardDeviation = standardDeviation*standardDeviation # sd^2
              standardDeviation = standardDeviation / (i) # i = (n-1)
              standardDeviation = math.sqrt(standardDeviation)
            else:
              standardDeviation = 100 # unrealisitic value to make sure we only calculate SD's for elements in distribution > 1

          # if the next point in the sample distribution is greater than cofactor * the running standard deviation,
          # accept it as the green-up day

          if((i+1) < len(halfSeason.intensities)):
            if( halfSeason.intensities[i+1] > (runningAverage + (cofactor * standardDeviation))):
              return halfSeason.days[i+1]

      # if we've considered the entire distribution and haven't found a spike in intensity, try some recursion and see if we can't pull out
      # a reasonable value by decreasing the cofactor.  If we can't decrease the cofactor anymore, return -1 as an error

      if(cofactor > 0):
        cofactor = cofactor-0.1
      else:
        return -1

  #
  # goldmanStaticValueMethod()
  # The Goldman static value method constructs a hard raster-intensity value out of training data.  It looks at the intensity
  # values of every sample distribution the user provided, and identifies the intensity value for each distribution at p=0.02.
  # It averages those intensity values and then sequentially searches each sample distribution for that value +/- FUZZ FACTOR,
  # returning the first day it finds of that value.

  def goldmanStaticValueMethod():
    return 0

  #
  # linearFit()
  # Preform a simple linear fit for all non-sense responseDist values (-1) that the pettorelliRunningAvgMethod
  # method failed to identify in the original distribution search.  This is experimental and assumes a linear
  # relationship between Petorelli and Goldman green-up/senescence days.  It is a heuristic tool that I slapped together
  # and is mentioned nowhere in the scientific literature (yet).

  def linearFit(self, y1, y2):
    
    # assuming a simple linear relationship of y1 = m(y2) + 0

    y1_scrubbed = list()
    y2_scrubbed = list()

    for var in y1:
      if var != -1:
        y1_scrubbed.append(var)

    for var in y2:
      if var != -1:
        y2_scrubbed.append(var)

    m = (sum(y1_scrubbed)/len(y1_scrubbed)) / (sum(y2_scrubbed)/len(y2_scrubbed))

