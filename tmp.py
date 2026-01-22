import paftacular as pft

print(pft.parse("y12{PEPTIDE-[Oxidation]}/3.4ppm*0.85")[0].mass(calculate_sequence=True))
print(pft.parse("y12{PEPTIDE-[Oxidation]}/3.4ppm*0.85")[0].dict_composition(calculate_sequence=True))
print(pft.parse("y12{P[+10]EPTIDE-[Oxidation]}/3.4ppm*0.85")[0].mass(calculate_sequence=True))
print(pft.parse("p-2[iTRAQ115]"))
print(pft.parse_single("p-2[iTRAQ115]+2[iTRAQ115][M-2H2O+Na]^2").as_dict())
print(pft.parse_single("p-2[iTRAQ115]+2[iTRAQ115][M-2H2O+Na]^2").ion_type)
