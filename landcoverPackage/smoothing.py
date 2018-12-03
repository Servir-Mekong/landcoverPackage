from utils import *

class Smoothing():

	def __init__(self):
		pass
	
	# add min max to the image
	def scale(self,img):
	  studyArea = img.geometry().bounds()
	  maxs = img.reduceRegion(reducer=ee.Reducer.percentile([99.99]),geometry=studyArea,scale=300,maxPixels=1e13)
	  mins = img.reduceRegion(reducer=ee.Reducer.percentile([0.01]),geometry=studyArea,scale=300,maxPixels=1e13)

	  return img.set("max",ee.Number(maxs.get("Mode"))).set("min",ee.Number(mins.get("Mode")))

	# scale the image
	def scaleImg(self,img):
	  maxs = ee.Number(img.get("max"))
	  mins = ee.Number(img.get("min"))
	  result = ee.Image(img.subtract(mins)).divide(maxs.subtract(mins)).multiply(100)
	 
	  return result.min(100).set('time_start',img.get('system:time_start'));
  

	# Function to compute the inverse log ratio of a regression results to
	# transform back to percent units
	def inverseLogRatio(self,image):
	  bands = image.bandNames()
	  ilrImage = ee.Image(100).divide(ee.Image(1).add(image.exp())).rename(bands)
	  return ilrImage

	def whittakerSmoothing(self,imageCollection, isCompositional = False, lamb = 5):
	 
		# quick configs to set defaults
		def toFl(image):
			return image.toFloat()

		# procedure start
		ic = imageCollection.map(toFl)

		dimension = ic.size()
		identity_mat = ee.Array.identity(dimension)
		difference_mat = getDifferenceMatrix(identity_mat,3)
		difference_mat_transpose = difference_mat.transpose()
		lamda_difference_mat = difference_mat_transpose.multiply(lamb)
		res_mat = lamda_difference_mat.matrixMultiply(difference_mat)
		hat_matrix = res_mat.add(identity_mat)

		# backing up original data
		original = ic

		def getProperties(image):
			return ee.Image(image).toDictionary()

		# get original image properties
		properties = ic.toList(10000).map(getProperties)


		# if data is compositional
		# calculate the logratio of an image between 0 and 100. First
		# clamps between delta and 100-delta, where delta is a small positive value.
		if (isCompositional):

			def clampImage(image):
				delta = 0.001
				bands = image.bandNames()
				image = image.clamp(delta,100-delta)
				image = (ee.Image.constant(100).subtract(image)).divide(image).log().rename(bands)
				return image

			ic = ic.map(clampImage)

		arrayImage = original.toArray()
		coeffimage = ee.Image(hat_matrix)
		smoothImage = coeffimage.matrixSolve(arrayImage)

		def getImageId(image):
			return ee.Image(image).id()

		idlist = ic.toList(10000).map(getImageId)

		bandlist = ee.Image(ic.first()).bandNames()

		flatImage = smoothImage.arrayFlatten([idlist,bandlist])
		smoothCollection = ee.ImageCollection(unpack(flatImage, idlist, bandlist))

		if (isCompositional):
			smoothCollection = smoothCollection.map(inverseLogRatio)

		def addSuffix(band):
			return ee.String(band).cat('_fitted')

		# get new band names by adding suffix fitted
		newBandNames = bandlist.map(addSuffix)

		# rename the bands in smoothened images
		smoothCollection = smoothCollection.select(bandlist, newBandNames)

		# a really dumb way to loose the google earth engine generated ID so that the two
		# images can be combined for the chart
		dumbimg = arrayImage.arrayFlatten([idlist,bandlist])
		dumbcoll = ee.ImageCollection(unpack(dumbimg,idlist, bandlist))
		outCollection = dumbcoll.combine(smoothCollection)

		outCollList = outCollection.toList(10000)
		def addPropBack(image):
			return ee.Image(image).set(properties.get(outCollList.indexOf(image)))

		outCollectionProp = outCollList.map(addPropBack)

		residue_sq = smoothImage.subtract(arrayImage).pow(ee.Image(2)).divide(dimension)
		rmse_array = residue_sq.arrayReduce(ee.Reducer.sum(),[0]).pow(ee.Image(1/2))

		rmseImage = rmse_array.arrayFlatten([["rmse"],bandlist])

		return (ee.ImageCollection(outCollectionProp), rmseImage)


	def setTime(self,img):
		return img.set("system:time_start",img.get("time_start"))



def whittakerSmoothen(collection):

	smooth = Smoothing()
	collection = collection.map(smooth.scale)
	collection = collection.map(smooth.scaleImg)
	smoothingResults, rmse = smooth.whittakerSmoothing(collection);
	smoothingResults = smoothingResults.map(smooth.setTime)

	return ee.ImageCollection(smoothingResults.select(["Mode_fitted"])), ee.Image(rmse)


	
	
