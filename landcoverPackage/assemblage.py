import ee

class assemblage():

	def __init__(self):
		pass
		
	def createAssemblage(self,image,nodeStruct):

		classStruct = self._toClassStruct(nodeStruct)
			
		# The starting id, i.e. the first decision
		startId = 'key1';
		
		# The initial decision tree string (DO NOT CHANGE)
		DTstring = ['1) root 9999 9999 9999'];
		# Call the function to construct the decision tree string (DO NOT CHANGE)
		DTstring = "\n".join(self.decision(nodeStruct,classStruct,startId,1,DTstring))#.join("\n");
		
		classifier = ee.Classifier.decisionTree(DTstring)
		
		sd = 10
		nIter = 100
		nBands = ee.List.sequence(0,image.bandNames().length().getInfo()-1,1)
		
	
		def monteCarlo(i):
			def createRand(j):
				rand = ee.Image(image.select(0)).gt(-1).multiply(ee.Image.random(ee.Number(i).multiply(j))).subtract(0.5).multiply(2)
				random = rand.multiply(sd)
				return ee.Image(image.select(ee.Number(j))).add(random)
			img = ee.ImageCollection(nBands.map(createRand))
			img = self.collectionToImage(img)
			classified = img.classify(classifier);
			return classified
		
		
		iters = ee.List.sequence(1,nIter)
		assemblage = ee.ImageCollection(iters.map(monteCarlo))
				
		mode = ee.Image(assemblage.mode())
		
		def uncertainty(img):
			return img.eq(mode)
		
		prob = ee.Image(assemblage.map(uncertainty).sum()).rename('prob')
				
		return mode, prob

	# Function to convert a dictionary of nodes into a decision tree string
	def decision(self,nodeStruct,classStruct,id1,node,DTstring):
		# Extract parameters
		lnode = 2*node; #// left child node number
		rnode = lnode + 1; #// right child node number
		dict1 = nodeStruct[id1]; #// current dictionary
		band = dict1['band']; #// decision band
		threshold = dict1['threshold']; #// decision threshold
		left = dict1['left']; #// left result (either 'terminal' or new id)
		right = dict1['right']; #// right result (either 'terminal' or new id)
		leftName = dict1['leftName']; #// left class name (if 'terminal')
		if right == 'terminal':
			rightName = dict1['rightName']; #// right class name (if 'terminal')
		
		leftLine = ''
		rightLine = '';
		leftNumber = 0;
		rightNumber = 0;

		# Add the left condition and right condition strings to the current decision 
		# tree string. If either condition is non-terminal, recursively call the function.
		if (left == 'terminal'): # left terminal condition
			leftNumber = str(classStruct[leftName]['number']);
			leftLine = str(lnode) + ') ' + str(band) + '>=' + str(threshold) + ' 9999 9999 ' + str(leftNumber) + ' *';
			DTstring.append(leftLine);
			if (right == 'terminal'): #// right terminal condition
				rightNumber = classStruct[rightName]['number'];
				rightLine = str(rnode) + ') ' + str(band) + '<' + str(threshold) + ' 9999 9999 ' + str(rightNumber) + ' *';
				DTstring.append(rightLine);
				return DTstring;
			else: # right non-terminal condition
				rightLine = str(rnode) + ') ' + str(band) + '<' + str(threshold) + ' 9999 9999 9999';
				DTstring.append(rightLine);
				return self.decision(nodeStruct,classStruct,right,rnode,DTstring);
    
		else: # left non-terminal condition
			leftLine = str(lnode) + ') ' + str(band) + '>=' + str(threshold) + ' 9999 9999 9999';
			DTstring.append(leftLine);	
			DTstring = self.decision(nodeStruct,classStruct,left,lnode,DTstring);
		
			if (right == 'terminal'): # // right terminal condition
				rightNumber = classStruct[rightName]['number'];
				rightLine = str(rnode) + ') ' + str(band) + '<' + str(threshold) + ' 9999 9999 ' + str(rightNumber) + ' *';
				DTstring.append(rightLine);
				return DTstring;
			else:  # right non-terminal
				rightLine = str(rnode) + ') ' + str(band) + '<' + str(threshold) + ' 9999 9999 9999';
				DTstring.append(rightLine);
				return self.decision(nodeStruct,classStruct,right,rnode,DTstring);
    
  
		return DTstring;

	def collectionToImage(self,collection):
		
		def iterate(img,prev):
			return ee.Image(prev).addBands(img)
		
		stack = ee.Image(collection.iterate(iterate,ee.Image(1)))
		 
		stack = stack.select(ee.List.sequence(1, stack.bandNames().size().subtract(1)));
		
		return stack;

	def _toClassStruct(self,nodeStruct):
		names = set([])
		for key, node in nodeStruct.iteritems():
			left_name = node.get('leftName')
			if left_name:
				names.add(left_name)
			right_name = node.get('rightName')
			if right_name:
				names.add(right_name)
		names.remove('other')
		names = list(names)
		names.sort()
		names.insert(0, 'other')

		class_struct = {}
		
		for i, name in enumerate(names):
			class_struct[name] = {'number': i}
		return class_struct


	

def createAssemblage(image,nodeStruct):
	
	assembl, uncertainty = assemblage().createAssemblage(image,nodeStruct)
	
	return assembl, uncertainty
	
