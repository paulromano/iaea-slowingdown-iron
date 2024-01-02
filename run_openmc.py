from pathlib import Path
import tarfile

import endf
import numpy as np
import openmc


cross_sections = [
    ('endfb80', '/lcrc/project/OpenMCValidation/data/hdf5/endfb-viii.0-hdf5/cross_sections.xml'),
    ('jeff33', '/lcrc/project/OpenMCValidation/data/hdf5/jeff-3.3-hdf5/cross_sections.xml'),
    ('tendl2021', '/lcrc/project/OpenMCValidation/data/hdf5/tendl-2021-hdf5/cross_sections.xml'),
    ('jendl5', '/lcrc/project/OpenMCValidation/data/hdf5/jendl-5.0-hdf5/cross_sections_e8photon.xml'),
]

# Iron composition: 100% Fe56, split by isotopes
# Source energy: 14 MeV, 2 MeV

photon_groups = 1e6*np.fromstring(Path('gamma-121gpr.txt').read_text(), sep=' ')
neutron_groups = 1e6*np.fromstring(Path('neutron-366gpr.txt').read_text(), sep=' ')

# Define materials
iron = openmc.Material(name='iron')
iron.add_nuclide('Fe56', 1.0)
iron.set_density('g/cm3', 7.8)

# Define geometry
inner_sphere1 = openmc.Sphere(r=10.0)
inner_sphere2 = openmc.Sphere(r=20.0)
boundary_sphere = openmc.Sphere(r=30.0, boundary_type='vacuum')
cell1 = openmc.Cell(fill=iron, region=-inner_sphere1)
cell2 = openmc.Cell(fill=iron, region=+inner_sphere1 & -inner_sphere2)
cell3 = openmc.Cell(fill=iron, region=+inner_sphere2 & -boundary_sphere)
model = openmc.Model()
model.geometry = openmc.Geometry([cell1, cell2, cell3])

# Define source
source_energy = 2.0e6
model.settings.source = openmc.IndependentSource(
    space=openmc.stats.Point(),
    energy=openmc.stats.Discrete([source_energy], [1.0])
)

# Define other settings
model.settings.run_mode = 'fixed source'
model.settings.particles = 10_000_000
model.settings.batches = 100
model.settings.photon_transport = True


# Define tallies
cell_filter = openmc.CellFilter([cell1, cell2, cell3])
neutron_filter = openmc.ParticleFilter(['neutron'])
neutron_energy = openmc.EnergyFilter(neutron_groups)
neutron_flux_tally = openmc.Tally(name='neutron flux')
neutron_flux_tally.filters = [neutron_filter, cell_filter, neutron_energy]
neutron_flux_tally.scores = ['flux']

photon_filter = openmc.ParticleFilter(['photon'])
photon_energy = openmc.EnergyFilter(photon_groups)
photon_flux_tally = openmc.Tally(name='photon flux')
photon_flux_tally.filters = [photon_filter, cell_filter, photon_energy]
photon_flux_tally.scores = ['flux']

model.tallies = [neutron_flux_tally, photon_flux_tally]

# Extract IRDFF-II ENDF file if it doesn't exist
if not Path('IRDFF-II.endf').is_file():
    with tarfile.open('IRDFF-II.tar.xz', 'r') as tar:
        print('Extracing IRDFF-II.tar.xz...')
        tar.extractall()

# Get tabulated cross sections from IRDFF-II
reactions = [
    ('Rh103', (10, 4)),
    ('In115', (10, 4)),
    ('Al27', (3, 107)),
    ('S32', (3, 103)),
]
endf_mats = endf.get_materials('IRDFF-II.endf', encoding='cp1250')
xs = {}
for nuclide, (mf, mt) in reactions:
    Z, A, m = endf.data.zam(nuclide)
    for mat in endf_mats:
        if mat[1, 451]['ZA'] == 1000*Z + A:
            if mf == 10:
                xs[nuclide] = mat[mf, mt]['levels'][0]['sigma']
            else:
                xs[nuclide] = mat[mf, mt]['sigma']

# Create spectral index tallies
for nuclide, sigma in xs.items():
    si_tally = openmc.Tally(name=f'Spectral index {nuclide}')
    multiplier = openmc.EnergyFunctionFilter(sigma.x, sigma.y)
    si_tally.filters = [neutron_filter, multiplier, cell_filter, neutron_energy]
    si_tally.scores = ['flux']
    model.tallies.append(si_tally)

# Create heating tally
heating_tally = openmc.Tally(name='heating')
all_particles = openmc.ParticleFilter(['neutron', 'photon', 'electron', 'positron'])
heating_tally.filters = [all_particles, cell_filter, neutron_energy]
heating_tally.scores = ['heating']
model.tallies.append(heating_tally)

# Create current tallies
neutron_current = openmc.Tally(name='neutron current')
surface_filter = openmc.SurfaceFilter([inner_sphere1, inner_sphere2, boundary_sphere])
neutron_current.filters = [neutron_filter, surface_filter, neutron_energy]
neutron_current.scores = ['current']
model.tallies.append(neutron_current)

photon_current = openmc.Tally(name='photon current')
photon_current.filters = [photon_filter, surface_filter, photon_energy]
photon_current.scores = ['current']
model.tallies.append(photon_current)

run_kwargs = {
    'mpi_args': ['mpiexec']
}

for label, xs_file in cross_sections:
    # Set cross sections
    openmc.config['cross_sections'] = xs_file

    # Case 1: Fe56, E=2.0 MeV
    sp_path = model.run(cwd=f'fe56_2MeV_{label}', **run_kwargs)

    # Case 2: Fe56, E=14.0 MeV
    model.settings.source[0].energy.x[0] = 14.0e6
    sp_path = model.run(cwd=f'fe56_14MeV_{label}', **run_kwargs)

    # Case 3: Fe, E=14.0 MeV
    iron.nuclides.clear()
    iron.add_nuclide('Fe54', 0.05845)
    iron.add_nuclide('Fe56', 0.91754)
    iron.add_nuclide('Fe57', 0.02119)
    iron.add_nuclide('Fe58', 0.00282)
    sp_path = model.run(cwd=f'fe_14MeV_{label}', **run_kwargs)

    # Case 4: Fe, E=2.0 MeV
    model.settings.source[0].energy.x[0] = 2.0e6
    sp_path = model.run(cwd=f'fe_2MeV_{label}', **run_kwargs)

