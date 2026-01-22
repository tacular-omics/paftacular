import paftacular as pf

print(pf.mzPAFParser().parse("y12/3.4ppm*0.85"))
print(pf.mzPAFParser().parse("p-2[iTRAQ115]"))
print(pf.mzPAFParser().parse_single("p-2[iTRAQ115]+2[iTRAQ115][M-2H2O+Na]^2").as_dict())
print(pf.mzPAFParser().parse_single("p-2[iTRAQ115]+2[iTRAQ115][M-2H2O+Na]^2").ion_type)
"""
{'ion_type': {}, 'analyte_reference': None, 'is_auxiliary': False, 'neutral_losses': 
[{'sign': -1, 'count': 2, '_formula': None, '_mass': None, '_reference': 'iTRAQ115'}, 
{'sign': 1, 'count': 2, '_formula': None, '_mass': None, '_reference': 'iTRAQ115'}], 
'isotopes': [], 'adducts': [{'sign': -1, 'count': 2, '_formula': 'H2O'}, {'sign': 1, 'count': 1, '_formula': 'Na'}], 
'charge': 2, 'mass_error': None, 'confidence': None}
"""
