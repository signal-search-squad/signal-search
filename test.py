#!/usr/bin/python

import sys
from itertools import chain, combinations

sys.path.append("models")

import sm
from sm import object_library as o

model = sm

def findByName(name, l):
	for p in l:
		if p.name == name:
			return p

def findParticleByName(name):
	return findByName(name, model.all_particles)

def findParameterByName(name):
	return findByName(name, model.all_parameters)

def toParticleList(l):
	l2 = []
	for c in l:
		if type(c) is str:
			l2.append(findParticleByName(c))
		else:
			l2.append(c)
	return l2


class ParticleSet:
	def __init__(self, l):
		self.l = l
	
	def __neg__(self):
		l2 = []
		for pi in self.l:
			if pi.name == pi.antiname:
				l2.append(pi)
			else:
				l2.append(findParticleByName(pi.antiname))
		return ParticleSet(l2)

class CompositeParticle(o.Particle):
	def __init__(self, name, components):
		o.Particle.__init__(self, pdg_code=999999,
			name = name,
			antiname = name + "~",
			spin = -1,
			color = -1,
			mass = model.parameters.ZERO,
			width = model.parameters.ZERO,
			texname = name,
			antitexname = name + "~",
			charge = 0,
			GhostNumber = 0,
			LeptonNumber = 0,
			Y = 0)
		self.components = ParticleSet(toParticleList(components))

p = CompositeParticle("p", ["g", "u", "u~", "d", "d~", "s", "s~", "c", "c~", "b", "b~"])
model.all_particles.append(p)
model.all_particles.append(p.anti())

vertices = {}

# make a dictionary to quickly find vertices containing a 
# certain particle
for v in model.all_vertices:
	for p in v.particles:
		if p in vertices:
			vertices[p].add(v)
		else:
			vertices[p] = set([v])

def findVerticesWithOutgoing(l):
	s = vertices[l[0]]
	for i in range(1, len(l)):
		s = s & vertices[l[i]]
	return s


def prettyPrintVertex(v, products):
	other = set(v.particles) - set(products)
	# at least one incoming particle
	# so we need all sub-sets of other with at least one element
	for i in range(1, len(other) + 1):
		for ss in combinations(other, i):
			ss = set(ss)
			s = ""
			for p in ss:
				s = s + p.name + " "
			s = s + " -> "
			rest = other - ss
			for p in rest:
				s = s + p.name + " "
			for p in products:
				s = s + p.name + " "
			print s

def diagrams(sources, products):
	sources = toParticleList(sources)
	products = toParticleList(products)
	v1 = findVerticesWithOutgoing(products)
	for v in v1:
		prettyPrintVertex(v, products)

diagrams(["p", "p~"], ["t", "t~"])