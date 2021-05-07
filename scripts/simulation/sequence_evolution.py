#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import random
import math


class sequenceEvol:
	'''
	Each instance of this class represens one evolutionary trajectory run. Starting with a sequence set of equal sequences,
	growing and mutating over time.
	'''

	def __init__(self, length, p_mut, N_init, t_start, t_final, p_repl=None, p_repl2=None, rand_repl=None, p_death=None, t_switch=None, init_seq=None):
		self.L = length
		self.p_repl = p_repl
		self.p_repl2 = p_repl2
		self.rand_repl = rand_repl
		self.p_death = p_death
		self.p_mut = p_mut
		self.N = N_init
		self.t_start = t_start
		self.t_final = t_final
		self.t_switch = t_switch
		self.init_seq = self._initSeq(init_seq)

	def _initSeq(self, init_seq=None):
		'''
		Create initial sequence as array of each nucleotide. If no sequence is given, create a random one.
		'''
		if init_seq is None:
			alphabet = 'ACGT'
			init_seq = ''.join(random.SystemRandom().choice(alphabet)for _ in range(self.L))
		return (init_seq)


	def evolve_poi(self):
		'''
		Mutate initial sequence set over the course of t generations
		'''
		print("--- Run simulation ---")

		# start time
		t = self.t_start+1

		# the current replication rate of the time step
		curr_p_repl = self.p_repl

		curr_N = self.N

		#list of dictionaries per time point, containing each sequence species an their amount
		species_dict_per_t = []
		species_dict_per_t.append({self.init_seq: curr_N})

		new_set = [self.init_seq]
		while t < self.t_final:
			print("Current time t = ", t)

			if self.rand_repl or self.p_repl is None:
				curr_p_repl = math.sin(t * 0.11) / 10 + 1.03
				#curr_p_repl = math.sin(t * 0.1) / 20 + 1.05
			else:
				# Switch the replication rate somewhere in the middle if that is the correct mode (p_repl2!=0)
				if (self.p_repl2 is not None) and (t >= self.t_switch):
					curr_p_repl = self.p_repl2

			print("Current replication rate = ", curr_p_repl)
			# draw number of sequences of the next generation
			curr_N = np.random.poisson(curr_p_repl * curr_N)

			print("Number of new sequences: ", curr_N)

			# collect sequence species for t
			species_dict = {}

			if curr_N > 0:
				new_set = random.choices(new_set, k=curr_N)

				# draw number of mutations in new generation
				# if number of new mutation sizes is larger than sites (rare/not likely) take number of sites
				# (= 2 mutations on one site, ignore back mutation)
				num_mut_sites = min(np.random.poisson(self.p_mut * curr_N * self.L), curr_N * self.L)
				print("Number of mutating sites: ",num_mut_sites)

				if num_mut_sites > 0:
					#for i in mut_sites:
					for _ in range(num_mut_sites):
						i = random.randrange(curr_N * self.L)
						seq_ind = math.floor(i / self.L)
						seq_pos = i % self.L
						seq = list(new_set[seq_ind])
						# seq[seq_pos] = self._mutateBase(new_set[seq_ind][seq_pos],initSeq[seq_pos])
						seq[seq_pos] = self.mutateBase(new_set[seq_ind][seq_pos])
						s = "".join(seq)
						new_set[seq_ind] = s

				# count species
				for s in new_set:
					species_dict[s] = species_dict.get(s, 0) + 1

			# add the sequences to the map
			species_dict_per_t.append(species_dict)
			t += 1

		return species_dict_per_t

	@staticmethod
	def mutateBase(base):
		'''
		Mutate individual nucleotide according
		to Kimura model, fixed transition probability
		'''

		stay = 0
		alpha = 0.2
		beta = 0.2
		gamma = 0.6

		transitions = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
		transversions_alpha = {'A': 'C', 'C': 'A', 'T': 'G', 'G': 'T'}
		transversions_beta = {'A': 'G', 'G': 'A', 'T': 'C', 'C': 'T'}

		r = np.random.uniform(0, 1)

		if stay <= r < (stay + gamma):
			return transitions[base]
		elif (stay + gamma) <= r < (stay + gamma + alpha):
			return transversions_alpha[base]
		elif (stay + gamma + alpha) <= r < (stay + gamma + alpha + beta):
			return transversions_beta[base]
		return base

