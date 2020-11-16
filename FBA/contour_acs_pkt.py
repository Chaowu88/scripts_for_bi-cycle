#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '06/02/2020'
__version__ = '1.0'




import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FormatStrFormatter
import cobra
from cobra import Model, Reaction, Metabolite




### setting
out_dir = r'path\to\output'
auto_model_file = r'path\to\iHN637_a.mat'
hetero_model_file = r'path\to\iHN637_h.mat'

auto_pkt_bnds = [0, 10]
auto_acs_bnds = [0, 10]
hetero_pkt_bnds = [0, 50]
hetero_acs_bnds = [0, 30]

n = 20


### simulation
os.makedirs(out_dir, exist_ok = True)

for model_file, mode, pkt_bnds, acs_bnds in zip([auto_model_file, hetero_model_file, hetero_model_file], 
												 ['Autotrophic mode', 'Heterotrophic mode', 'Mixotrophic mode'], 
												 [auto_pkt_bnds, hetero_pkt_bnds, hetero_pkt_bnds], 
												 [auto_acs_bnds, hetero_acs_bnds, hetero_acs_bnds]): 
	
	# build model
	model = cobra.io.load_matlab_model(model_file)
	
	xu5p = model.metabolites.get_by_id('xu5p-D[c]')
	pi = model.metabolites.get_by_id('pi[c]')
	gap = model.metabolites.get_by_id('g3p[c]')
	actp = model.metabolites.get_by_id('actp[c]')

	pkt = Reaction('pkt')
	pkt.add_metabolites({xu5p: -1, pi: -1, gap: 1, actp: 1})
	model.add_reactions([pkt])
	
	
	r_etoh = model.reactions.get_by_id('EX_etoh(e)')
	r_ac = model.reactions.get_by_id('EX_ac(e)')
	model.objective = {r_etoh: 1, r_ac: 1}
	
	
	if mode == 'Autotrophic mode':
		model.reactions.get_by_id('EX_h2(e)').lower_bound = -10
		model.reactions.get_by_id('EX_h2(e)').upper_bound = 0 
		model.reactions.get_by_id('EX_co(e)').lower_bound = -10
		model.reactions.get_by_id('EX_co(e)').upper_bound = 0 
		model.reactions.get_by_id('EX_co2(e)').lower_bound = -10
		model.reactions.get_by_id('EX_co2(e)').upper_bound = 0 

	elif mode == 'Heterotrophic mode':
		model.reactions.get_by_id('EX_fru(e)').lower_bound = -100
		model.reactions.get_by_id('EX_fru(e)').upper_bound = 0
		model.reactions.get_by_id('EX_co2(e)').lower_bound = 0
		model.reactions.get_by_id('EX_co2(e)').upper_bound = 100

	else:
		model.reactions.get_by_id('EX_co(e)').lower_bound = -100
		model.reactions.get_by_id('EX_co(e)').upper_bound = 0
		model.reactions.get_by_id('EX_fru(e)').lower_bound = -100
		model.reactions.get_by_id('EX_fru(e)').upper_bound = 0
	
	
	# simulate
	x = np.linspace(pkt_bnds[0], pkt_bnds[1], n)
	y = np.linspace(acs_bnds[0], acs_bnds[1], n)
	X, Y = np.meshgrid(x, y)

	Z = np.empty((n, n))
	for i in range(n):
		for j in range(n):
			
			with model:
				cons1 = model.problem.Constraint(model.reactions.get_by_id('pkt').flux_expression, lb = X[i, j], ub = X[i, j])   
				cons2 = model.problem.Constraint(model.reactions.get_by_id('CODH_ACS').flux_expression, lb = Y[i, j], ub = Y[i, j])  
				
				model.add_cons_vars(cons1)
				model.add_cons_vars(cons2)
				
				if mode != 'Autotrophic mode':
					cons3 = model.problem.Constraint(model.reactions.get_by_id('PFL').flux_expression, lb = 0, ub = 0)
					cons4 = model.problem.Constraint(model.reactions.get_by_id('HCO3E').flux_expression, lb = 0, ub = 0)
					
					model.add_cons_vars(cons3)
					model.add_cons_vars(cons4)
				
				model.optimize()
				
				Z[i, j] = model.reactions.get_by_id('EX_etoh(e)').flux + model.reactions.get_by_id('EX_ac(e)').flux
	Z[Z < 0] = 0
	
	
	# plot
	plt.figure()

	plt.contourf(X, Y, Z, 100, cmap = plt.cm.get_cmap('RdBu').reversed())

	ax = plt.gca()
	ax.tick_params(labelsize = 25)
	plt.xlabel('pkt flux (mmol gDCW$^{-1}$h$^{-1}$)', fontsize = 25)
	plt.ylabel('acs flux (mmol gDCW$^{-1}$h$^{-1}$)', fontsize = 25)

	cbar = plt.colorbar()
	cbarTicks = cbar.get_ticks()
	newCbarTicks = np.linspace(cbarTicks.min(), cbarTicks.max(), 5)
	if mode == 'Autotrophic mode':
		cbar.set_ticks(newCbarTicks.round(1))
	else:
		cbar.set_ticks(newCbarTicks.round(0))
	cbar.ax.tick_params(labelsize = 25)
	cbar.set_label('EtOH and AC production rate\n(mmol gDCW$^{-1}$h$^{-1}$)', labelpad = 60, rotation = 270, fontsize = 25)

	plt.suptitle(mode, x = 0.5, y = 1, fontsize = 30)

	plt.savefig('%s/%s.jpg' % (out_dir, mode), dpi = 300, bbox_inches = 'tight')
	plt.close()
	
	
	





