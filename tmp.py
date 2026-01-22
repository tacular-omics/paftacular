import paftacular as pft

print(pft.parse("y12{PEPTIDE}/3.4ppm*0.85")[0].mass(calculate_sequence=True))
print(pft.parse("p-2[iTRAQ115]"))
print(pft.parse_single("p-2[iTRAQ115]+2[iTRAQ115][M-2H2O+Na]^2").as_dict())
print(pft.parse_single("p-2[iTRAQ115]+2[iTRAQ115][M-2H2O+Na]^2").ion_type)
"""
{'ion_type': {}, 'analyte_reference': None, 'is_auxiliary': False, 'neutral_losses': 
[{'sign': -1, 'count': 2, '_formula': None, '_mass': None, '_reference': 'iTRAQ115'}, 
{'sign': 1, 'count': 2, '_formula': None, '_mass': None, '_reference': 'iTRAQ115'}], 
'isotopes': [], 'adducts': [{'sign': -1, 'count': 2, '_formula': 'H2O'}, {'sign': 1, 'count': 1, '_formula': 'Na'}], 
'charge': 2, 'mass_error': None, 'confidence': None}
"""
