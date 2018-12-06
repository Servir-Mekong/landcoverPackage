import ee, math
import landsat
import covariates
import smoothing
import assemblage


def composite(aoi,year,sensors="Landat8"):
	# PARAM AOI: Area of interest (Feature)
	# PARAM Year: (integer)
	# PARAM Sensors : (Landat4,Landat5,Landat7,Landsat8 )
	
	# function returns an image

	image = landsat.composite(aoi,year)
	return image


def addCovariates(img):
	# function returns an image

	image = covariates.returnCovariates(img)
	return image



def sample(composite,trainingData):
	# PARAM composite image with (blue, green, red, nir, swir1, thermal, swir2)  
	# PARAM training data FeatureCollection
	
	# RETURN FeatureCollection
	
	composite = addCovariates(composite)
	training = ee.FeatureCollection(composite.sampleRegions(trainingData, ["class"], 30))		
	 
	return training


def primitive(composite,primitiveType,trainingData,year):
	# PARAM composite image with (blue, green, red, nir, swir1, thermal, swir2)  
	# PARAM primitive type
	# PARAM training data FeatureCollection
	
	# RETURN image

	composite = addCovariates(composite)
	classifier = ee.Classifier.randomForest(100).setOutputMode('PROBABILITY').train(trainingData,"class")
	classification = composite.classify(classifier,'Mode')
	
	classification = classification.set({"year":year})
  
	return ee.Image(1).subtract(classification).multiply(100).toInt()


def smoothing(primitive):
	# PARAM primitives: imagecollection  
	
	# RETURN image collection and rmse

	collection, rmse = smoothing.whittakerSmoothen(primitive)
	
	return collection, rmse

def getAssemblage(image,nodeStruct):
	# PARAM primitives: image with primitives as bands  
	# Node structure see example
	#  	nodeStruct = { 	'key1':  {'band': 'aquaculture','threshold': 50, 'left': 'terminal', 'leftName': 'aquaculture', 'right': 'key2'},
	#					'key2':  {'band': 'barren', 'threshold': 40, 'left': 'terminal', 'leftName': 'barren', 'right': 'key3'},
	#					'key3':  {'band': 'cropland', 'threshold': 60, 'left': 'terminal', 'leftName': 'cropland', 'right': 'key4'},
	#					'key4':  {'band': 'forest', 'threshold': 5, 'left': 'terminal', 'leftName': 'other', 'right': 'terminal', 'rightName': 'forest'}	};
	#
	# RETURN assemblage: image: numbers represent classes in alphabetical order starting with 0 = other
	# RETURN probability: image: result from monte carlo
	
	assem, prob = assemblage.createAssemblage(image,nodeStruct)
	
	return assem, prob

