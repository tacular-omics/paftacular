import paftacular as pft

print(pft.parse_multi("y12{PEPTIDE-[Oxidation]}/3.4ppm*0.85")[0].mass(calculate_sequence=True))
print(pft.parse_multi("y12{PEPTIDE-[Oxidation]}/3.4ppm*0.85")[0].dict_composition(calculate_sequence=True))
print(pft.parse_multi("y12{P[+10]EPTIDE-[Oxidation]}/3.4ppm*0.85")[0].mass(calculate_sequence=True))
print(pft.parse_multi("p-2[iTRAQ115]"))
print(pft.parse_single("p-2[iTRAQ115]+2[iTRAQ115][M-2H2O+Na]^2").as_dict())
print(pft.parse_single("p-2[iTRAQ115]+2[iTRAQ115][M-2H2O+Na]^2").ion_type)


print(pft.PrecursorIon().__repr__())
ion1 = pft.ImmoniumIon(pft.AminoAcids.K, None)
ion2 = pft.ImmoniumIon(pft.AminoAcids.R, None)

print(ion1._cache is ion2._cache)  # True - same object
print(ion1._cache is pft.ImmoniumIon._cache)  # True - class variable
