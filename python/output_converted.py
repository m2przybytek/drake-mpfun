import glob
import os
import re


def convert_output(file_name=None):
    if not os.path.isdir("raw"):
        os.mkdir("raw")

    if not os.path.isdir("results"):
        os.mkdir("results")

    for file_ in [file_name] if file_name else glob.glob("*.out"):
        with open("raw/" + file_, "r") as infile:
            result = infile.read()

        scf_energy = re.findall("SCF energy:\s+([-\d.]+)", result)
        mp2_energy = re.findall("Total MP2 energy \(standard\)\s*\n\s*([-\d.]+)", result)
        mp2_ec_energy = re.findall("FINAL RESULTS : TOTAL MP2\s*\n.*?\n\s+.*?\s+([-\d.]+)", result)
        fccd_energy = re.findall("FINAL RESULTS : TOTAL FCCD\s*\n.*?\n\s+.*?\s+([-\d.]+)", result)

        if scf_energy:
            with open(f"results/{file_}_SCF", "w") as outfile:
                outfile.write(scf_energy[0])

        if mp2_energy:
            with open(f"results/{file_}_MP2", "w") as outfile:
                outfile.write(mp2_energy[0])

        if mp2_ec_energy:
            with open(f"results/{file_}_MP2_EC", "w") as outfile:
                outfile.write(mp2_ec_energy[0])

        if fccd_energy:
            with open(f"results/{file_}_FCCD", "w") as outfile:
                outfile.write(fccd_energy[0])


if __name__ == "__main__":
    convert_output()
