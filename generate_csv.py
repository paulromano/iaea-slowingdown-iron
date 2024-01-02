import openmc
import pandas as pd


libraries = ['endfb80', 'jeff33', 'jendl5', 'tendl2021']
source_energies = ['14MeV', '2MeV']
composition = ['fe', 'fe56']

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


def get_results(composition, source, library):
    with openmc.StatePoint(f'{composition}_{source}_{library}/statepoint.100.h5') as sp:
        for name in names:
            tallies[name] = sp.get_tally(name=name)

    neutron_energy_lower = tallies['neutron flux'].filters[-1].bins[:, 0]
    neutron_energy_upper = tallies['neutron flux'].filters[-1].bins[:, 1]
    neutron_flux = tallies['neutron flux'].get_reshaped_data().squeeze()
    neutron_current = tallies['neutron current'].get_reshaped_data().squeeze()
    spectral_rh103 = tallies['Spectral index Rh103'].get_reshaped_data().squeeze()
    spectral_in115 = tallies['Spectral index In115'].get_reshaped_data().squeeze()
    spectral_al27 = tallies['Spectral index Al27'].get_reshaped_data().squeeze()
    spectral_s32 = tallies['Spectral index S32'].get_reshaped_data().squeeze()
    heating = tallies['heating'].get_reshaped_data().squeeze().sum(axis=0)
    neutron_df = pd.DataFrame({
        'E_lo [eV]': neutron_energy_lower,
        'E_hi [eV]': neutron_energy_upper,
        'Neutron flux (0-10cm)': neutron_flux[0],
        'Neutron flux (10-20cm)': neutron_flux[1],
        'Neutron flux (20-30cm)': neutron_flux[2],
        'Neutron current (10cm)': neutron_current[0],
        'Neutron current (20cm)': neutron_current[1],
        'Neutron current (30cm)': neutron_current[2],
        'Rh103 spectra index (0-10cm)': spectral_rh103[0],
        'Rh103 spectra index (10-20cm)': spectral_rh103[1],
        'Rh103 spectra index (20-30cm)': spectral_rh103[2],
        'In115 spectra index (0-10cm)': spectral_in115[0],
        'In115 spectra index (10-20cm)': spectral_in115[1],
        'In115 spectra index (20-30cm)': spectral_in115[2],
        'Al27 spectra index (0-10cm)': spectral_al27[0],
        'Al27 spectra index (10-20cm)': spectral_al27[1],
        'Al27 spectra index (20-30cm)': spectral_al27[2],
        'S32 spectra index (0-10cm)': spectral_s32[0],
        'S32 spectra index (10-20cm)': spectral_s32[1],
        'S32 spectra index (20-30cm)': spectral_s32[2],
        'Total heating (0-10cm)': heating[0],
        'Total heating (10-20cm)': heating[1],
        'Total heating (20-30cm)': heating[2],
    })
    neutron_df.to_csv(f'{composition}_{source}_{library}_neutron.csv', index=False)

    photon_energy_lower = tallies['photon flux'].filters[-1].bins[:, 0]
    photon_energy_upper = tallies['photon flux'].filters[-1].bins[:, 1]
    photon_flux = tallies['photon flux'].get_reshaped_data().squeeze()
    photon_current = tallies['photon current'].get_reshaped_data().squeeze()
    photon_df = pd.DataFrame({
        'E_lo [eV]': photon_energy_lower,
        'E_hi [eV]': photon_energy_upper,
        'Photon flux (0-10cm)': photon_flux[0],
        'Photon flux (10-20cm)': photon_flux[1],
        'Photon flux (20-30cm)': photon_flux[2],
        'Photon current (10cm)': photon_current[0],
        'Photon current (20cm)': photon_current[1],
        'Photon current (30cm)': photon_current[2],
    })
    photon_df.to_csv(f'{composition}_{source}_{library}_photon.csv', index=False)


for comp in composition:
    for E in source_energies:
        for library in libraries:
            get_results(comp, E, library)

