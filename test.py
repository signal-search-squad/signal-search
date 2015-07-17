#!/usr/bin/python

import sys
from itertools import chain, combinations, product, permutations, izip

sys.path.append("models")

import sm as model
from sm import object_library as o

import UFOExtensions
UFOExtensions.patch(o, model)

from lib import *

import logging

logging.basicConfig(filename = "log.txt", format = "%(asctime)-15s %(message)s")
log = logging.getLogger("main")
log.setLevel(logging.INFO)


# thank you Mr. Knuth
def algorithm_u(ns, m):
	def visit(n, a):
		ps = [[] for i in xrange(m)]
		for j in xrange(n):
			ps[a[j + 1]].append(ns[j])
		return ps

	def f(mu, nu, sigma, n, a):
		if mu == 2:
			yield visit(n, a)
		else:
			for v in f(mu - 1, nu - 1, (mu + sigma) % 2, n, a):
				yield v
		if nu == mu + 1:
			a[mu] = mu - 1
			yield visit(n, a)
			while a[nu] > 0:
				a[nu] = a[nu] - 1
				yield visit(n, a)
		elif nu > mu + 1:
			if (mu + sigma) % 2 == 1:
				a[nu - 1] = mu - 1
			else:
				a[mu] = mu - 1
			if (a[nu] + sigma) % 2 == 1:
				for v in b(mu, nu - 1, 0, n, a):
					yield v
			else:
				for v in f(mu, nu - 1, 0, n, a):
					yield v
			while a[nu] > 0:
				a[nu] = a[nu] - 1
				if (a[nu] + sigma) % 2 == 1:
					for v in b(mu, nu - 1, 0, n, a):
						yield v
				else:
					for v in f(mu, nu - 1, 0, n, a):
						yield v

	def b(mu, nu, sigma, n, a):
		if nu == mu + 1:
			while a[nu] < mu - 1:
				yield visit(n, a)
				a[nu] = a[nu] + 1
			yield visit(n, a)
			a[mu] = 0
		elif nu > mu + 1:
			if (a[nu] + sigma) % 2 == 1:
				for v in f(mu, nu - 1, 0, n, a):
					yield v
			else:
				for v in b(mu, nu - 1, 0, n, a):
					yield v
			while a[nu] < mu - 1:
				a[nu] = a[nu] + 1
				if (a[nu] + sigma) % 2 == 1:
					for v in f(mu, nu - 1, 0, n, a):
						yield v
				else:
					for v in b(mu, nu - 1, 0, n, a):
						yield v
			if (mu + sigma) % 2 == 1:
				a[nu - 1] = 0
			else:
				a[mu] = 0
		if mu == 2:
			yield visit(n, a)
		else:
			for v in b(mu - 1, nu - 1, (mu + sigma) % 2, n, a):
				yield v

	n = len(ns)
	a = [0] * (n + 1)
	for j in xrange(1, m + 1):
		a[n - m + j] = j - 1
	return f(m, n, 0, n, a)


def subsets(n, set):
	if n == 1:
		return iter([[set]])
	else:
		return algorithm_u(set, n)


def findByName(name, l):
	for p in l:
		if p.name == name:
			return p


def findParticleByName(name):
	return findByName(name, model.all_particles)


def findParameterByName(name):
	return findByName(name, model.all_parameters)


def antiparticle(p):
	if p.antiname == p.name:
		return p
	else:
		return findParticleByName(p.antiname)


def antiparticles(l):
	return [antiparticle(p) for p in l]


def toParticleList(l):
	l2 = []
	for c in l:
		if type(c) is str:
			l2.append(findParticleByName(c))
		else:
			l2.append(c)
	return l2


def toInstanceList(l):
	return [ParticleInstance(p) for p in l]

def particleTypes(l):
	return [p.particleType for p in l]

# We need this because some of the UFO stuff uses negation
# on all the fields in a class to get the anti-particle,
# so when that is done on CompositeParticle below and it
# hits the sub-particle list, it fails if it's a plain list
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

	def anti(self):
		return self.__neg__()


class ParticleFractions:
	def __init__(self, l):
		self.l = l

	def __neg__(self):
		return self

	def __getitem__(self, item):
		return self.l[item]

class CompositeParticle(o.Particle):
	def __init__(self, name, components, fractions, antiname = None):
		o.Particle.__init__(self, pdg_code=999999,
			name = name,
			antiname = name + "~" if antiname == None else antiname,
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
		self.fractions = ParticleFractions(fractions)

	def anti(self):
		cp = CompositeParticle(self.antiname, [], self.fractions, antiname = self.name)
		cp.components = self.components.anti()
		return cp

	def getAmplitude(self, p):
		# TODO: simple PDFs
		for i in range(0, len(self.components.l)):
			if p == self.components.l[i]:
				return self.fractions.l[i]
		return 0.0

p = CompositeParticle("p",
	["g", "u", "u~", "d", "d~", "s", "s~",  "c", "c~", "b",    "b~"],
	[0.5, 0.6, 0.01, 0.3, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001])
model.all_particles.append(p)
model.all_particles.append(p.anti())

vertices = {}

# make a dictionary to quickly find vertices containing a 
# certain particle
for v in model.all_vertices:
	# skip ghosts
	# TODO: why?
	hasGhosts = False
	for p in v.particles:
		if p.GhostNumber != 0:
			hasGhosts = True
			break
	if hasGhosts:
		continue
	if len(set(v.particles)) == 1:
		# ignore things like ggg, gggg
		continue
	for p in v.particles:
		if p in vertices:
			vertices[p].add(v)
		else:
			vertices[p] = set([v])


# A vertex that knows the difference between outgoing and incoming
class OrderedVertex:
	def __init__(self, incoming, outgoing, coupling, additionalOutgoing):
		self.incoming = incoming
		self.outgoing = outgoing
		self.coupling = coupling
		self.additionalOutgoing = additionalOutgoing
	
	def __repr__(self):
		return str(self.incoming) + " -> " + str(self.outgoing) + "  (" + str(self.coupling) + ")"
	
	def isNull(self):
		return False


class NullVertex(OrderedVertex):
	def __init__(self, particleInstance):
		self.incoming = [particleInstance]
		self.outgoing = [particleInstance]
		self.coupling = 1.0
		self.additionalOutgoing = []
	
	def isNull(self):
		return True

def getAmplitude(p1, p2):
	if p1 == p2:
		return 1.0
	elif isinstance(p2, CompositeParticle):
		return p2.getAmplitude(p1)
	else:
		return 0.0

def listsEqual(l1, l2):
	return MultiSet(particleTypes(l1)) == MultiSet(particleTypes(l2))

class Diagram:
	def __init__(self, particles, amplitude, prev = None, vertices = None):
		self.particles = particles
		self.amplitude = amplitude
		self.prev = prev
		self.vertices = vertices
	
	def copy(self):
		nd = Diagram([], self.amplitude, self.prev, self.vertices)
		nd.particles = nd.particles + self.particles
		return nd

	def getOverlap(self, sources):
		if len(self.particles) == 1:
			return 0.0
		if len(sources) < len(self.particles):
			return 0.0
		# consider all possible permutations of sources
		# since, for example, a p p~ pair has two ways
		# of producing a u u~ pair (PDFs for both u and u~ are non-zero for both p and p~)
		amplitude = 0.0
		perms = permutations(sources)
		for perm in perms:
			permAmplitude = 1.0
			match = True
			for pair in izip(self.particles, perm):
				a = getAmplitude(pair[0].particleType, pair[1])
				if a == 0.0:
					match = False
					break
				else:
					permAmplitude = permAmplitude * a
			if match:
				log.info("Match for %s %s", self.particles, perm)
				amplitude = amplitude + permAmplitude
		return amplitude

	def hasLoops(self):
		if self.prev == None or self.prev.vertices == None:
			return False
		for v1 in self.vertices:
			for v2 in self.prev.vertices:
				if listsEqual(v2.outgoing, v1.incoming) and listsEqual(v1.outgoing, v2.incoming):
					return True
		return False

	def __repr__(self):
		return "Diargam[amplitude = " + str(self.amplitude) + self.__repr2__()

	def __repr2__(self):
		return "\n\t" + str(self.particles) + "   |   " + str(self.vertices) + ("" if self.prev == None else self.prev.__repr2__())

# a class that works like a set but allows duplicates that it keeps count of
class MultiSet:
	def __init__(self, iterable):
		self.d = {}
		self.len = 0
		for i in iterable:
			self.add(i)

	def add(self, x):
		if x in self.d:
			self.d[x] = self.d[x] + 1
		else:
			self.d[x] = 1
		self.len = self.len + 1

	def remove(self, x, count = 1):
		if x in self.d:
			n = self.d[x]
			if n <= count:
				del self.d[x]
				self.len = self.len - n
			else:
				self.d[x] = n - count
				self.len = self.len - count


	def __sub__(self, o):
		copy = MultiSet([])
		copy.d = self.d.copy()
		copy.len = self.len
		for key in o.d:
			copy.remove(key, count = o.d[key])
		return copy

	def __iter__(self):
		l = []
		for key in self.d:
			l = l + [key] * self.d[key]
		return l.__iter__()

	def __repr__(self):
		return self.d.__repr__()

	def __eq__(self, other):
		if self.len != other.len:
			return False
		return self.d == other.d

	def __len__(self):
		return self.len


def findVerticesWithOutgoing(l):
	types = particleTypes(l)
	sTypes = MultiSet(types)
	s = vertices[l[0].particleType]
	for i in range(1, len(l)):
		s = s & vertices[l[i].particleType]
	
	r = []
	# now go through all possible arrangements of outgoing/incoming
	for v in s:
		allParticles = MultiSet(v.particles)
		remaining = allParticles - sTypes

		if len(remaining) > 0:
			# one choice is that all remaining are incoming
			r.append(OrderedVertex(toInstanceList(antiparticles(remaining)), l, v.getCoupling(), []))

		# need at least one incoming, so generate all 2-partitions 
		part = subsets(2, list(remaining))
		for p in part:
			if len(p[0]) == 0:
				# skip the no-incoming particle results
				continue
			incoming2 = antiparticles(p[0])
			if MultiSet(incoming2) == sTypes:
				# if just propagating and radiating something, ignore
				continue
			extras = toInstanceList(p[1])
			r.append(OrderedVertex(toInstanceList(incoming2), extras + l, v.getCoupling(), extras))
	
	if len(l) == 1:
		# just propagating is an option
		r.append(NullVertex(l[0]))
	return r


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

def found(d):
	print "Found:"
	print d
	log.info("Found %s", d)

def filterLastStepVertices(depth, vs, maxSrc):
	if depth == 1:
		return [v for v in vs if len(v.incoming) <= maxSrc]
	else:
		return vs

def diagrams(sources, partialDiagram, amplitudeThreshold, depth):
	if depth == 0:
		log.info("Maxdepth reached")
		return
	intermediate = partialDiagram.particles
	
	for nsubsets in range(1, len(intermediate) + 1):
		partitions = subsets(nsubsets, intermediate)
		
		# There is some intermediate state that contains some particles
		# Each of them, and various combinations of them can come from
		# various vertices. So divide the intermediate particles into
		# all possible sub-sets. This is a partition. Then for each partition,
		# go through all diagrams that can produce each sub-set in the partition
		# and move to the next step
		
		for partition in partitions:
			log.info("Partition: %s", partition)
			# each subset is a set of intermediate particles
			# than need to be part of a single vertex
			# for each subset, find all vertices that could have
			# the subset as outgoing particles
			vertexChoices = []
			partLen = len(partition)
			# each vertex has at least one incoming
			# so the maximum number of incoming particles
			# that this vertex can have is the total minus one
			# for each other vertex
			maxVertexIncoming = len(sources) - partLen + 1
			for si in range(0, partLen):
				subset = partition[si]

				vertexChoices.append(filterLastStepVertices(depth, findVerticesWithOutgoing(subset), maxVertexIncoming))
			
			# pick all possible combinations of diagrams
			# for example, if we have the following intermediate particles:
			#
			# a, b, c
			#
			# and the following vertices match:
			#
			# a1, a2 -> a
			# a3, a4 -> a
			# b1, b2 -> b
			# b3, b4 -> b
			# c1, c2 -> c
			#
			# then we need to explore the following possibilities:
			# 
			# [a1, a2 -> a, b1, b2 -> b, c1, c2 -> c]
			# [a1, a2 -> a, b3, b4 -> b, c1, c2 -> c]
			# [a3, a4 -> a, b1, b2 -> b, c1, c2 -> c]
			# [a3, a4 -> a, b3, b4 -> b, c1, c2 -> c]
			#
			# so we generate 4 new sets of intermediate particles:
			# 
			# [a1, a2, b1, b2, c1, c2], [a1, a2, b3, b4, c1, c2], [a3, a4, b1, b2, c1, c2], [a3, a4, b3, b4, c1, c2]
			# 
			# and recursively explore each
			log.info("Current: %s", partialDiagram)
			log.info("Vertices:")
			for v in vertexChoices:
				log.info("\t[")
				for vv in v:
					log.info("\t\t%s", vv)
				log.info("\t]")
			
			for p in product(*vertexChoices):
				i2 = []
				v2 = []
				amplitude = partialDiagram.amplitude
				# copy to allow extending intermediate list without
				# affecting other branches
				pdcopy = partialDiagram.copy()
				allNull = True
				for v in list(p):
					amplitude = amplitude * v.coupling / 4.0
					i2 = i2 + v.incoming
					v2.append(v)
					pdcopy.particles = pdcopy.particles + v.additionalOutgoing
					if not listsEqual(v.incoming, v.outgoing) or len(v.additionalOutgoing) > 0:
						allNull = False

				if log.isEnabledFor(logging.DEBUG):
					log.debug("Considering %s", v2)
					log.debug("Current amplitude: %s", amplitude)
				if allNull:
					log.debug("Discarding (just propagating) %s", v2)
					continue
				
				# check if below amplitude threshold
				if abs(amplitude) < amplitudeThreshold:
					log.debug("Discarding (low amplitude) %s", v2)
					continue
				newpd = Diagram(i2, amplitude, pdcopy, v2)

				# check for loop diagrams
				if newpd.hasLoops():
					log.debug("Discarding (loop) %s", v2)
					continue


				a = newpd.getOverlap(sources)
				if a > 0.0:
					found(newpd)
				log.info("Recursing: ")
				log.info("\t%s", pdcopy.particles)
				log.info("\t%s", newpd.particles)
				diagrams(sources, newpd, amplitudeThreshold, depth - 1)

				
		#print partition
		#v1 = findVerticesWithOutgoing(intermediate)
		
		#for v in v1:
		#	prettyPrintVertex(v, products)

incoming = toParticleList(["p", "p~"])
outgoing = toInstanceList(toParticleList(["t", "t~"]))

cd = Diagram(outgoing, 1.0)
diagrams(incoming, cd, 0.001, 2)
#diagrams([], [1, 2, 3], 0.1)
