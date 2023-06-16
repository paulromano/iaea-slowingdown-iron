import openmc
import argparse
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('statepoint')
args = parser.parse_args()

nuclides = ['Rh103', 'In115', 'Al27', 'S32']

tallies = {}
names = [
    'neutron flux',
    'photon flux',
    'Spectral index Rh103',
    'Spectral index In115',
    'Spectral index Al27',
    'Spectral index S32',
    'heating',
    'neutron current',
    'photon current'
]

with openmc.StatePoint(args.statepoint) as sp:
    for name in names:
        tallies[name] = sp.get_tally(name=name)

neutron_energy = tallies['neutron flux'].filters[1].values
neutron_flux = tallies['neutron flux'].mean.ravel()

photon_energy = tallies['photon flux'].filters[1].values
photon_flux = tallies['photon flux'].mean.ravel()

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
ax1.stairs(neutron_flux, neutron_energy, label='Neutron')
ax1.stairs(photon_flux, photon_energy, label='Photon')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('Energy [eV]')
ax1.set_ylabel('Flux')
ax1.legend()

# Plot spectral indices
for nuc in nuclides:
    spectral_index = tallies[f'Spectral index {nuc}'].mean.ravel()
    ax2.stairs(spectral_index, neutron_energy, label=nuc)
ax2.set_yscale('log')
ax2.set_xlim(0., 15.e6)
ax2.set_xlabel('Energy [eV]')
ax2.set_ylabel('Spectral index')
ax2.legend(loc='lower right')

# Get heating as 2D array with shape (4, groups)
heating = tallies['heating'].get_reshaped_data().squeeze()
neutron_heating = heating[0]
photon_heating = heating[1:].sum(axis=0)
total_heating = heating.sum(axis=0)

ax3.stairs(neutron_heating, neutron_energy, label='Neutron')
ax3.stairs(photon_heating, neutron_energy, label='Photon')
ax3.stairs(total_heating, neutron_energy, label='Total')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel('Energy [eV]')
ax3.set_ylabel('Heating')
ax3.legend()

# Plot neutron and photon current at 10, 20, and 30 cm
neutron_current = tallies['neutron current'].get_reshaped_data().squeeze()
ax4.stairs(neutron_current[0], neutron_energy, label='Neutron, 10cm')
ax4.stairs(neutron_current[1], neutron_energy, label='Neutron, 20cm')
ax4.stairs(neutron_current[2], neutron_energy, label='Neutron, 30cm')
photon_current = tallies['photon current'].get_reshaped_data().squeeze()
ax4.stairs(photon_current[0], photon_energy, label='Photon, 10cm')
ax4.stairs(photon_current[1], photon_energy, label='Photon, 20cm')
ax4.stairs(photon_current[2], photon_energy, label='Photon, 30cm')
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xlabel('Energy [eV]')
ax4.set_ylabel('Current')
ax4.set_ylim(ymin=1e-6)
ax4.legend()
plt.show()
