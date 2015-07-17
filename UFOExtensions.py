#!/usr/bin/python

import sys
import importlib

def particleRepr(o):
	return o.name

# Also patch to get a numeric coupling

def getCoupling2(o):
	# this is wrong, but this is it for now
	if not hasattr(o, "_numericCoupling"):
		expr = o.couplings[(0, 0)].value
		try:
			o._numericCoupling = eval(expr)
		except:
			print expr
			raise
	return o._numericCoupling

class EvalDict:
	def __init__(self, _dict):
		self._dict = _dict
		self.values = {}
	
	def __len__(self):
		pass
	
	def __getitem__(self, key):
		if key in self.values:
			return self.values[key]
		elif key in self._dict:
			obj = self._dict[key]
			if type(obj).__name__ == "Parameter":
				# if numeric, cache and return it
				# if a string expression, recursively evaluate it
				value = obj.value
				if type(value).__name__ == "str":
					value = eval(value, {}, self)
					self.values[key] = value
					print key + " = " + str(value)
					return value
				else:
					self.values[key] = value
					return value
			else:
				self.values[key] = obj
				return obj
		else:
			raise KeyError("Variable not found: '" + key + "'")
			
	
	def __setitem__(self, key, value):
		pass
	
	def __delitem__(self, key):
		pass
	
	def __iter__(self):
		pass
	
	def __contains__(self, key):
		pass

def patch(o, model):
	print "Patching " + str(o)
	sys.path.append("models/" + model.__name__)
	
	params = __import__("parameters", globals(), locals(), [], -1)
	functs = __import__("function_library", globals(), locals(), [], -1)
	
	_dict = params.__dict__.copy()
	_dict.update(functs.__dict__)
	
	ed = EvalDict(_dict)
	
	
	def vertex_getCoupling(o):
		# this is wrong, but this is it for now
		if not hasattr(o, "_numericCoupling"):
			expr = o.couplings[(0, 0)].value
			try:
				o._numericCoupling = eval(expr, {}, ed)
			except:
				print expr
				raise
		return o._numericCoupling
	
	def parameter_getNumericValue(p):
		if not hasattr(o, "_numericValue"):
			expr = o.value
			try:
				o._numericValue = eval(expr, {}, ed)
			except:
				print expr
				raise
	
	# when printing a Particle object, make it print its name
	o.Particle.__repr__ = particleRepr
	o.Vertex.getCoupling = vertex_getCoupling
	o.Parameter.getNumericValue = parameter_getNumericValue
	