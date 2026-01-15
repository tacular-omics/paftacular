import paftacular as pf

print(pf.mzPAFParser().parse("y12/3.4ppm*0.85"))
print(pf.mzPAFParser().parse("p-2[iTRAQ115]"))
print(pf.mzPAFParser().parse("p-2[iTRAQ115]+2[iTRAQ115][M-2H2O+Na]^2"))
