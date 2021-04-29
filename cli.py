import pydicom
import sys
import logging
import argparse

import indexcalc

parser = argparse.ArgumentParser(description="Calcula indices de complejidad para un plan de VMAT/IMRT")
parser.add_argument("infile", type=str, help="Fichero de entrada con el RTPlan")
parser.add_argument("--debug", action="store_true", default=False,
                    help="Activa generación de fichero debug.txt"
)
args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG, filename="debug.txt", filemode="w")
else:
    logging.basicConfig(level=logging.ERROR)

fname = args.infile
data = pydicom.read_file(fname, force=True)

mu_per_gy = indexcalc.mu_per_gy(data)
print(f"MU/Gy = {mu_per_gy:.2f}")

sas = indexcalc.sas(data, 2)
print(f"SAS(2mm) = {sas:.2f}%")

plan_pi = indexcalc.pi(data)
print(f"PI = {plan_pi:.2f}")