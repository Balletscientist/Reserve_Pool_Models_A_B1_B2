#!/usr/bin/env python3

from sre_constants import ANY_ALL
import sys
import os
import math

MCELL_PATH = os.environ.get('MCELL_PATH', '')
if MCELL_PATH:
     sys.path.append(os.path.join(MCELL_PATH, 'lib'))
else:
     print("Error: variable MCELL_PATH that is used to find the mcell library was not set.")
     sys.exit(1)

import mcell as m

ITERATIONS = 1000

cyt =  m.geometry_utils.create_icosphere(
    name = 'CYT', 
    radius = 0.35, 
    subdivisions = 4
)

# ---- observables ----
SEED = 1
viz_output = m.VizOutput(
     mode = m.VizMode.ASCII,
     output_files_prefix = './viz_data/seed_' + str(SEED).zfill(5) + '/Scene',
     every_n_timesteps = 1
)

observables = m.Observables()
observables.add_viz_output(viz_output)

model = m.Model()

model.add_geometry_object(cyt)

CaMKII_species = m.Species(
    name = 'CaMKII',
    diffusion_constant_3d = 1e-6
)
model.add_species(CaMKII_species)

Actin_species = m.Species(
    name = 'Actin',
    diffusion_constant_3d = 0.000000001
)
model.add_species(Actin_species)

Mobile_species = m.Species(
    name = 'Mobile',
    diffusion_constant_3d = 5e-5
)
model.add_species(Mobile_species)


rel_a = m.ReleaseSite(
     name = 'rel_a',
     complex = Actin_species,
     region = cyt,
     number_to_release = 500
)
model.add_release_site(rel_a)

rel_c = m.ReleaseSite(
     name = 'rel_c',
     complex = CaMKII_species,
     shape = m.Shape.SPHERICAL,
     location = (0, -0.34, 0),
     number_to_release = 1000
)
model.add_release_site(rel_c)

count_mobile = m.Count (
    expression = m.CountTerm(
        species_pattern = Mobile_species,
    ),
    file_name = './mobilised_SV_no_exclusion.dat',
)

observables.add_count(count_mobile)

model.add_observables(observables)

# Binding reaction rate - mobilising synaptic vesciles
k_onCaMKII = 5e7

# Callback function invoked every time there is a binding event - optional
def rxn_callback(rxn_info, context):
    context.count += 1
    print("Bound: " + str(context.count));

# create reaction rule
rxn_rule = m.ReactionRule(
    name = "actin binding",
    reactants=[ Actin_species, CaMKII_species ], 
    products=[ Mobile_species, CaMKII_species ],
    fwd_rate = k_onCaMKII,
)
model.add_reaction_rule(rxn_rule)
print(rxn_rule.to_bngl_str())

# thats for the call back to see the message - when the reaction happens - optional 
class RxnCallbackContext():
    def __init__(self):
        self.count = 0

context = RxnCallbackContext()

model.config.total_iterations = ITERATIONS
model.config.time_step = 1e-5
# model.config.seed = 12

model.initialize()
# model.dump_internal_state()

# Add the callback to log every time a binding event happens - optional 
model.register_reaction_callback(
    rxn_callback, context, rxn_rule 
)

for frame_number in range(1, ITERATIONS + 1):
    model.export_viz_data_model()
            
    model.run_iterations(1)


model.end_simulation()

# $MCELL_PATH/utils/visualize.sh viz_data/seed_00001/ 