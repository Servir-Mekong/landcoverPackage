import ee


def composite():
	## create the composite
	image = ee.Image("LANDSAT/LC08/C01/T1_SR/LC08_133048_20140113")
	image = image.select(ee.List([1,2,3,4,5,7,6,9,10,11],ee.List(['blue','green','red','nir','swir1','thermal','swir2'])))
	print image.bandNames.getInfo()

def training():
	## create the trainin gdata
	print "say hello"



def primitive():
	## create the primitive
	print "bye"



ee.Initialize()
#composite()
