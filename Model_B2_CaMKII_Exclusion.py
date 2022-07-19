#!/usr/bin/env python3

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

pool =  m.geometry_utils.create_icosphere(
    name = 'POOL', 
    radius = 0.3,
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
model.add_geometry_object(pool)

CaMKII_species = m.Species(
    name = 'CaMKII',
    diffusion_constant_3d = 1e-6
)
model.add_species(CaMKII_species)

Actin_species = m.Species(
    name = 'Actin',
    diffusion_constant_3d = 0.000000001 #3.1e-12
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
     region = pool,
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
    file_name = './camkii_exclusion.dat',
)

observables.add_count(count_mobile)

model.add_observables(observables)

# Binding reaction rate
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
model.config.seed = 23

transp = m.SurfaceClass(
     name = 'transp',
     type = m.SurfacePropertyType.TRANSPARENT,
     affected_complex_pattern = Actin_species
)
model.add_surface_class(transp)
pool.surface_class = transp

model.initialize()
# model.dump_internal_state()

# Add the callback to log every time a binding event happens - optional 
model.register_reaction_callback(
    rxn_callback, context, rxn_rule 
)
#and that will shrink by 3 quaters from the initial size. the end frame is quater of the size of the start frame. 
shrink_amount = ITERATIONS*4/3

for frame_number in range(1, ITERATIONS + 1):
    model.export_viz_data_model()
            
    model.run_iterations(1)
    
    expand_factor = 1 / (shrink_amount + 1 - frame_number)
    
    for vertex_number in range(len(pool.vertex_list)):
      vertex = pool.vertex_list[vertex_number]
      model.add_vertex_move(pool, vertex_number, (
          - vertex[0] * expand_factor,
          - vertex[1] * expand_factor,
          - vertex[2] * expand_factor ))

    model.apply_vertex_moves()

# for frame_number in range(1, ITERATIONS + 1):
#     model.export_viz_data_model()
            
#     model.run_iterations(1)

#     expand_factor = (frame_number + 1) / (frame_number * 500)
    
#     for vertex_number in range(len(pool.vertex_list)):
#       vertex = pool.vertex_list[vertex_number]
#       model.add_vertex_move(pool, vertex_number, (
#           - vertex[0] * expand_factor,
#           - vertex[1] * expand_factor,
#           - vertex[2] * expand_factor ))

#     model.apply_vertex_moves()

model.end_simulation()

# $MCELL_PATH/utils/visualize.sh viz_data/seed_00001/