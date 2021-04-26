import pydicom
import sys

import indexcalc


fname = sys.argv[1]
data = pydicom.read_file(fname, force=True)
mu_per_gy = indexcalc.mu_per_gy(data)
print(f"MU/Gy = {mu_per_gy:.2f}")

sas = indexcalc.sas(data, 2)
print(f"SAS(2mm) = {sas:.2f}%")