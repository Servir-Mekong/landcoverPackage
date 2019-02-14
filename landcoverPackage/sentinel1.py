import ee
import math 
from utils import *

class env(object):
	
	def __init__(self):
		"""Initialize the environment"""
		self.studyArea = ee.FeatureCollection("projects/Sacha/AncillaryData/StudyRegions/Ecuador_EcoRegions_Complete")
		self.polar = ["VV","VH"]
		self.direction = 'ASCENDING'; 
		self.startDate = ee.Date.fromYMD(2017,1,1)
		self.endDate = ee.Date.fromYMD(2017,12,31)
		
		self.percentiles = [20,80]
		
		self.applyMaskLowEntropy = True
		self.speckleFilter = True
		self.calculateRatio = True
	
		self.harmonics = 1
		
		self.ksize = 3
		self.enl = 7; 
		

class functions():       
	def __init__(self):
		"""Initialize the Surfrace Reflectance app."""  
	
	    # get the environment
		self.env = env() 	

	def main(self):
		collection = self.getSAR()
				
		if  self.env.applyMaskLowEntropy == True:
			collection = collection.map(self.maskLowEntropy)
			
				
		if self.env.speckleFilter == True:
			collection = collection.map(self.applySpeckleFilter)
			
			
		if self.env.calculateRatio == True and len(self.env.polar) >1:
			collection = collection.map(self.addRatio)
			
		
		
		perc = collection.reduce(ee.Reducer.percentile(self.env.percentiles))
		stddev = collection.reduce(ee.Reducer.stdDev())		
		
		# fix the harmonics here
		#self.getHarmonics(collection,'VV')
		
		
		
	def getSAR(self):
		""" get the Sentinel 1 data collection"""
		# Inventory of S1 dual polarization images.
		collection = ee.ImageCollection('COPERNICUS/S1_GRD').filterBounds(self.env.studyArea)\
														   .filter(ee.Filter.eq('instrumentMode','IW'))\
														   .filter(ee.Filter.eq('orbitProperties_pass', self.env.direction))\
														   .filter(ee.Filter.eq('transmitterReceiverPolarisation', self.env.polar))\
														   .filter(ee.Filter.date(self.env.startDate, self.env.endDate));
		

		
		
		
		return collection


	def maskLowEntropy(self,image): 
		""" Removes low-entropy edges """
		bad = image.select(0).multiply(10000).toInt().entropy(ee.Kernel.circle(5)).lt(3.2)
		return image.updateMask(image.mask().multiply(bad.focal_max(5).Not()))

	def applySpeckleFilter(self,img):
		vv = img.select('VV')
		vh = img.select('VH')
		vv = self.speckleFilter(vv).rename('VV');
		vh = self.speckleFilter(vh).rename('VH');
		return ee.Image.cat(vv,vh).copyProperties(img,['system:time_start']);


	def speckleFilter(self,image):
		""" apply the speckle filter """
		# Convert image from dB to natural values
		nat_img = self.toNatural(image);

		# Square kernel, ksize should be odd (typically 3, 5 or 7)
		weights = ee.List.repeat(ee.List.repeat(1,self.env.ksize),self.env.ksize);

		# ~~(ksize/2) does integer division in JavaScript
		kernel = ee.Kernel.fixed(self.env.ksize,self.env.ksize, weights, ~~(self.env.ksize/2), ~~(self.env.ksize/2), False);

		# Get mean and variance
		mean = nat_img.reduceNeighborhood(ee.Reducer.mean(), kernel);
		variance = nat_img.reduceNeighborhood(ee.Reducer.variance(), kernel);

		# "Pure speckle" threshold
		ci = variance.sqrt().divide(mean);# square root of inverse of enl

		# If ci <= cu, the kernel lies in a "pure speckle" area -> return simple mean
		cu = 1.0/math.sqrt(self.env.enl);

		# If cu < ci < cmax the kernel lies in the low textured speckle area
		# -> return the filtered value
		cmax = math.sqrt(2.0) * cu;

		alpha = ee.Image(1.0 + cu*cu).divide(ci.multiply(ci).subtract(cu*cu));
		b = alpha.subtract(self.env.enl + 1.0);
		d = mean.multiply(mean).multiply(b).multiply(b).add(alpha.multiply(mean).multiply(nat_img).multiply(4.0*self.env.enl));
		f = b.multiply(mean).add(d.sqrt()).divide(alpha.multiply(2.0));

		# If ci > cmax do not filter at all (i.e. we don't do anything, other then masking)

		# Compose a 3 band image with the mean filtered "pure speckle", 
		# the "low textured" filtered and the unfiltered portions
		out = ee.Image.cat(self.toDB(mean.updateMask(ci.lte(cu))),self.toDB(f.updateMask(ci.gt(cu)).updateMask(ci.lt(cmax))),image.updateMask(ci.gte(cmax)));	
		
		return out.reduce(ee.Reducer.sum());

	def addRatio(self,img):
		vv = self.toNatural(img.select(['VV'])).rename(['VV']);
		vh = self.toNatural(img.select(['VH'])).rename(['VH']);
		ratio = vh.divide(vv).rename(['ratio']);
		return ee.Image.cat(vv,vh,ratio).copyProperties(img,['system:time_start']);


	def toNatural(self,img):
		"""Function to convert from dB to natural"""
		return ee.Image(10.0).pow(img.select(0).divide(10.0));
		
	def toDB(self,img):
		""" Function to convert from natural to dB """
		return ee.Image(img).log10().multiply(10.0);


	def getHarmonics(self,collection,dependent): #function getHarmonicsAmplitudePhase(collection,dependent,harmonics,roi){
		#######################################/
		# Function that performs a harmonic time series analysis of an image 
		# collection for a dependent variable and returns a multi-band image of the  
		# harmonic coefficients and RMSE.
		
		# Make a list of harmonic frequencies to model.  
		# These also serve as band name suffixes.
		harmonicFrequencies = range(self.env.harmonics);
		harmonicSeq = ee.List.sequence(1,self.env.harmonics);
		  
		# Function to get a sequence of band names for harmonic terms.
		
		def getNames(base,myList):
			names = []
			for i in myList:
				names.append(base+str(i+1))
			return names
			#def getList(i):
			
			#	return ee.String(base).cat(ee.String(i))
			#return myList.map(getList)
			
		# Construct lists of names for the harmonic terms.
		cosNames = getNames('cos_', harmonicFrequencies);
		sinNames = getNames('sin_', harmonicFrequencies);
		amplitudeNames = getNames('amplitude_', harmonicFrequencies);
		phaseNames = getNames('phase_', harmonicFrequencies);


		# Independent variables.
		independents = ee.List(['constant', 't']).cat(cosNames).cat(sinNames);
		
		
		def addConstant(image):
			return image.addBands(ee.Image(1));

		# Function to add a time band.
		def addTime(image):
			# Compute time in fractional years since the epoch.
			date = ee.Date(image.get('system:time_start'));
			years = date.difference(ee.Date('1970-01-01'), 'year');
			timeRadians = ee.Image(years.multiply(2 * math.pi));
			return image.addBands(timeRadians.rename('t').float());
		  
		  
		# Function to compute the specified number of harmonics	
		# and add them as bands.  Assumes the time band is present.
		def addHarmonics(image):
			
			def getImg(freqs):
				# Make an image of frequencies.
				frequencies = ee.Image.constant(freqs);
				# This band should represent time in radians.
				time = ee.Image(image).select('t');
				# Get the cosine terms.
				cosines = time.multiply(frequencies).cos().rename(cosNames);
				# Get the sin terms.
				sines = time.multiply(frequencies).sin().rename(sinNames);
				return ee.Image(image.addBands(cosines).addBands(sines))
			
			for i in harmonicFrequencies:
				image = getImg(i)
			
			return image
			
		  
		# Add constants, time, and harmonics to image
		harmonicImages = collection.map(addConstant).map(addTime).map(addHarmonics);
		
		# The output of the regression reduction is a 4x1 array image.
		harmonicTrend = harmonicImages.select(independents.add(dependent)).reduce(ee.Reducer.linearRegression(independents.length(), 1));

		# Turn the array image into a multi-band image of coefficients.
		harmonicTrendCoefficients = harmonicTrend.select('coefficients').arrayProject([0]).arrayFlatten([independents]);
		harmonicTrendResiduals = harmonicTrend.select('residuals').arrayGet(0).rename('RMSE');	  

		out = harmonicTrendCoefficients.select(['constant','t']);
		  
		# Compute phase and amplitude
		def getPhaseAmp(harmonic):
			harmonicString = ee.String(ee.Number(harmonic).int());
			cosName_i = ee.String('cos_').cat(harmonicString);
			sinName_i = ee.String('sin_').cat(harmonicString);
			phaseName_i = ee.String('phase_').cat(harmonicString);
			amplitudeName_i = ee.String('amplitude_').cat(harmonicString);
			phase = harmonicTrendCoefficients.select(cosName_i).atan2(harmonicTrendCoefficients.select(sinName_i));
			phase = phase.rename(phaseName_i);
			
			amplitude = harmonicTrendCoefficients.select(cosName_i).hypot(harmonicTrendCoefficients.select(sinName_i));
			amplitude = amplitude.rename(amplitudeName_i);
			out = phase.addBands(amplitude);
			return out;
		
		phaseAmplitude = harmonicSeq.map(getPhaseAmp)
				
		out = out.addBands(self.newCollectionToImage(phaseAmplitude));
		
		print(out.bandNames().getInfo())
		out = out.addBands(harmonicTrendResiduals);

		# Rename band to dependent_band
		bands = ee.Image(out).bandNames();
		
		print(bands.getInfo())
		
	def newCollectionToImage(self,collection):

		""" Helper function to convert image collection into stack of image bands """
		def stacking(img,prev):	
			 return ee.Image(prev).addBands(img);
			
		stack = ee.Image(collection.iterate(stacking,ee.Image(1)))
		stack = stack.select(ee.List.sequence(1, stack.bandNames().size().subtract(1)));
	
		return stack


		
		
if __name__ == "__main__":
	ee.Initialize()
	functions().main()
