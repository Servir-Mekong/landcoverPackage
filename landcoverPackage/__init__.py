import ee, math
import landsat
import covariates


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
	training = ee.FeatureCollection(composite.sampleRegions(trainingData, ["land_class"], 30))		
	 
	return training


def primitive(composite,primitiveType,trainingData,year):
	# PARAM composite image with (blue, green, red, nir, swir1, thermal, swir2)  
	# PARAM primitive type
	# PARAM training data FeatureCollection
	
	# RETURN image

	classifier = ee.Classifier.randomForest(100).setOutputMode('PROBABILITY').train(trainingData,"class")
	classification = composite.classify(classifier,'Mode')
	
	classification = classification.set({"year":year})
  
	return ee.Image(1).subtract(classification).multiply(100).toInt()


def assemblage(year,primitiveCollection,):
	print "I'm just a dummy"
	img = ee.Image(primitiveCollection.first())
	
	return img
