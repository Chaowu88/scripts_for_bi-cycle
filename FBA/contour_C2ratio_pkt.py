#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'
__date__ = '06/01/2020'
__version__ = '1.0'




import os
import numpy as np
import matplotlib.pyplot as plt
import cobra
from cobra import Model, Reaction, Metabolite




### setting
out_dir = r'path\to\output'
model_file = r'path\to\iHN637_h.mat'

pkt_bnds = [0, 20]
ratio_bnds = [0.2, 2]   # etoh/ac

n = 20


### simulation
os.makedirs(out_dir, exist_ok = True)

# build model
model = cobra.io.load_matlab_model(model_file)

xu5p = model.metabolites.get_by_id('xu5p-D[c]')
pi = model.metabolites.get_by_id('pi[c]')
gap = model.metabolites.get_by_id('g3p[c]')
actp = model.metabolites.get_by_id('actp[c]')

pkt = Reaction('pkt')
pkt.add_metabolites({xu5p: -1, pi: -1, gap: 1, actp: 1})
model.add_reactions([pkt])

model.reactions.get_by_id('CODH_ACS').knock_out()

# simulate
x = np.linspace(pkt_bnds[0], pkt_bnds[1], n)
y = np.linspace(ratio_bnds[0], ratio_bnds[1], n)
X, Y = np.meshgrid(x, y)

Z = np.empty((n, n))
for i in range(n):
	for j in range(n):
		
		with model:
			cons1 = model.problem.Constraint(model.reactions.get_by_id('EX_etoh(e)').flux_expression - Y[i, j]*model.reactions.get_by_id('EX_ac(e)').flux_expression, lb = 0, ub = 0)
			cons2 = model.problem.Constraint(model.reactions.get_by_id('pkt').flux_expression, lb = X[i, j], ub = X[i, j]) 
			cons3 = model.problem.Constraint(model.reactions.get_by_id('PFL').flux_expression, lb = 0, ub = 0)
			
			model.add_cons_vars(cons1)
			model.add_cons_vars(cons2)
			model.add_cons_vars(cons3)
			
			model.optimize()
			
			Z[i, j] = model.reactions.get_by_id('EX_etoh(e)').flux + model.reactions.get_by_id('EX_ac(e)').flux

Z[Z < 0] = 0

# plot
plt.figure()

plt.contourf(X, Y, Z, 100, cmap = plt.cm.get_cmap('RdBu').reversed())

ax = plt.gca()
ax.tick_params(labelsize = 25)
plt.xlabel('pkt flux (mmol gDCW$^{-1}$h$^{-1}$h)', fontsize = 25)
plt.ylabel('Production flux Ratio of\nEtOH to AC', fontsize = 25)

cbar = plt.colorbar()
cbar.ax.tick_params(labelsize = 25)
cbarTicks = cbar.get_ticks()
newCbarTicks = np.linspace(cbarTicks.min(), cbarTicks.max(), 5)
cbar.set_ticks(newCbarTicks.round(1))
cbar.set_label('EtOH and AC production rate\n(mmol gDCW$^{-1}$h$^{-1}$)', labelpad = 60, rotation = 270, fontsize = 25)

plt.savefig('%s/c2ratio_pkt.jpg' % out_dir, dpi = 300, bbox_inches = 'tight')
plt.close()





