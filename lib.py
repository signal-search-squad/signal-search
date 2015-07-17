#!/usr/bin/python

class ParticleInstance:
	sid = 1
	def __init__(self, particleType):
		self.particleType = particleType
		ParticleInstance.sid = ParticleInstance.sid + 1
		self.id = ParticleInstance.sid
	
	def __repr__(self):
		return str(self.particleType) + "_" + str(self.id)
	
