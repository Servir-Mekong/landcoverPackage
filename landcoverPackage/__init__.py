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

	training = ee.FeatureCollection(composite.sampleRegions(trainingData, ["land_class"], 30))		
	 
	return training


def primitive(composite,primitiveType,trainingData,year):
	# PARAM composite image with (blue, green, red, nir, swir1, thermal, swir2)  
	# PARAM primitive type
	# PARAM training data FeatureCollection
	
	# RETURN image

	classifier = ee.Classifier.randomForest(100).setOutputMode('PROBABILITY').train(trainingData,"land_class")
	classification = composite.classify(classifier,'Mode')
	
	classification = classification.set({"year":year})
  
	return ee.Image(1).subtract(classification).multiply(100).toInt()



if __name__ == "__main__":  
	
	ee.Initialize()
	aoi = ee.Geometry.Polygon([[103.876,18.552],[105.806,18.552],[105.806,19.999],[103.876,19.999],[103.876,18.552]])
	trainingData = ee.FeatureCollection("users/servirmekong/NorthVietnam/urban2") #.randomColumn("rand")
	
	trainingData = trainingData #.filter(ee.Filter.gt("rand",0.50))
	print trainingData.size().getInfo()
	img = ee.Image("users/servirmekong/temp/0luse009")
	
	year = 2014
	
	#img = ee.Image(composite(aoi,year))
	#print "---------------"
	#print img.bandNames().getInfo()
	
	img = addCovariates(img)
	#print "---------------"
	#print img.bandNames().getInfo()
	
	samples = ee.FeatureCollection(sample(img,trainingData))
	#print samples.first().getInfo()
	primi = primitive(img,"rice",samples,2014)

	#print primi.bandNames().getInfo()

	task_ordered= ee.batch.Export.image.toAsset(image=primi, 
								  description="test", 
								  assetId="users/servirmekong/temp/0luse019" ,
								  region=aoi['coordinates'], 
								  maxPixels=1e13,
								  scale=30)	

	task_ordered.start() 
	
	#def getComposite(primitiveCollection,year):
		# PARAM primitive collection
		# RETURN assemblage
		# RETURN probability
		
	#	primi1 = ee.Image(primitiveCollection.first())
	#	primi2 = primi1.add(ee.Image.random())
	
	#	return primi1, primi2

	#def getCovariates

