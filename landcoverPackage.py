import ee


def composite(aoi,year,sensor="Landat"):
	# function takes in Area of interest (Feature)
	# function takes in Year (integer)
	# function takes in name of the Sensor
	
	# function returns an image
	
	image = ee.Image("LANDSAT/LC08/C01/T1_SR/LC08_133048_20140113")
	image = image.select(ee.List([1,2,3,4,5,7,8]),ee.List(['blue','green','red','nir','swir1','thermal','swir2']))
	
	return image

def training():
	## create the trainin gdata
	print "say hello"



def primitive():
	## create the primitive
	print "bye"



ee.Initialize()
composite()
