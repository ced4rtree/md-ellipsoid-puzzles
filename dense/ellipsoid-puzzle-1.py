#!/usr/bin/env python
# coding: utf-8

# In[52]:


import warnings

warnings.filterwarnings("ignore")


# In[53]:


from flowermd.library import EllipsoidChain, EllipsoidForcefield
from flowermd import Pack, Simulation
from flowermd.utils import create_rigid_ellipsoid_chain, get_target_box_number_density
import unyt as u
from cmeutils.gsd_utils import ellipsoid_gsd
import gsd.hoomd


# In[54]:


ellipsoids = EllipsoidChain(lengths=1, num_mols=128, bead_mass=1.0, lpar=1.0)


# In[55]:


system = Pack(
    density=0.15*u.Unit("nm**-3"),
    molecules=ellipsoids,
    packing_expand_factor=6,
    edge=2, overlap=1,
    fix_orientation=True
)


# In[56]:


forcefield = EllipsoidForcefield(
    bond_k=100, bond_r0=0,
    epsilon=1.0,
    lpar=1.0,
    lperp=0.5,
    r_cut=2.0
)
rigid_frame, rigid_constraint = create_rigid_ellipsoid_chain(snapshot=system.hoomd_snapshot)


# In[57]:


sim = Simulation(
    constraint=rigid_constraint,
    forcefield=forcefield.hoomd_forces,
    gsd_file_name="ellipsoid_density.gsd",
    log_file_name="ellipsoid_density.txt",
    initial_state=rigid_frame,
    gsd_write_freq=100
)


# Initial Parameters:
#  - density: 1.0
#  - period: 10
#  - n_steps: 6e4
#  - kT: 7.0
#  - tau_kt: 100*sim.dt
#  - thermalize_particles: True
# 
# Parameter changelog:
#  - density 5.0 -> fail
#  - density 4.0 -> fail
#  - density 3.0 -> fail
#  - kT 5.0 (from 7.0) -> fail
#  - kT 1.0 -> fail
#  - kT 10.0 -> fail
#  - density 2.0
#  - kT 0.0001 -> fail
#  - kT 5.0 & n_steps 6e5 -> fail
#  - n_steps 6e4 & period 100 -> fail
#  - period 1 -> fail
#  - period 100000 -> success! ...but nothing happens in the gsd?
#  - n_steps 5e6

# In[58]:


target_box = get_target_box_number_density(2.0*u.Unit("nm**-3"), n_beads=128)
sim.run_update_volume(
    final_box_lengths=target_box,
    period=100000,
    n_steps=5e6,
    kT=5.0,
    tau_kt=100*sim.dt,
    thermalize_particles=True,
)


# In[41]:


sim.run_NVT(n_steps=2e4, kT=1.0, tau_kt=0.001)
sim.flush_writers()


# In[44]:


ellipsoid_gsd_better(gsd_file="ellipsoid_density.gsd", new_file="ovito-ellipsoid_density.gsd", ellipsoid_types='R', lpar=1.0, lperp=0.5)


# In[43]:


def ellipsoid_gsd_better(gsd_file, new_file, ellipsoid_types, lpar, lperp):
    """Add needed information to GSD file to visualize ellipsoids.

    Saves a new GSD file with lpar and lperp values populated
    for each particle. Ovito can be used to visualize the new GSD file.

    Parameters
    ----------
    gsd_file : str
        Path to the original GSD file containing trajectory information
    new_file : str
        Path and filename of the new GSD file
    ellipsoid_types : str or list of str
        The particle types (i.e. names) of particles to be drawn
        as ellipsoids.
    lpar : float
        Value of lpar of the ellipsoids
    lperp : float
        Value of lperp of the ellipsoids

    """
    with gsd.hoomd.open(new_file, "w") as new_t:
        with gsd.hoomd.open(gsd_file) as old_t:
            for snap in old_t:
                shape_dicts_list = []
                for ptype in snap.particles.types:
                    if ptype == ellipsoid_types or ptype in ellipsoid_types:
                        shapes_dict = {
                            "type": "Ellipsoid",
                            "a": lpar,
                            "b": lperp,
                            "c": lperp,
                        }
                    else:
                        shapes_dict = {"type": "Sphere", "diameter": 0.001}
                    shape_dicts_list.append(shapes_dict)
                snap.particles.type_shapes = shape_dicts_list
                snap.validate()
                new_t.append(snap)


# In[ ]:




