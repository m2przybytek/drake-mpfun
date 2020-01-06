import datetime
import os
import subprocess
import time

import output_converted

import numpy as np
import termcolor


def two_digits(float_number):
    return "%.2f" % float_number


def construct_file_name(nocc, nbas, omega, gfac, alpha, c0, c1, bold=False):
    """Before was f'nocc={nocc}nbas={nbas}omega={omega}gfac={gfac}alpha={alpha}c0={c0}c1={c1}.out'"""

    def part(name, value):
        if bold:
            return f"{name}=" + termcolor.colored(f"{value}", attrs=["bold"])
        return f"{name}={value}"

    printout = (
        part("nocc", nocc)
        + part("nbas", nbas)
        + part("omega", omega)
        + part("gfac", gfac)
        + part("alpha", alpha)
        + part("c0", c0)
        + part("c1", c1)
    )

    if bold:
        return printout
    return printout + ".out"


def current_time():
    return datetime.datetime.now().isoformat(" ", timespec="seconds")


example_input_file = {
    "calc_type": "SCF",
    "nocc": "2",
    "nbas": "6",
    "omega": "6",  # oba indeksy niezależnie od zera do omegi phi_i phi_j
    "gfac": "2.",
    "alpha": "1./4.",
    "c0": "1.",  # c0 * 1
    "c1": "1./2.",  # c1 * |x|
    # tylko stosunek pomiędzy c0 i c1 jest istotnym
    "mandatory_keyword": "END",
    "max_and_min_eta": "4 6",  # wymusza silną ortogonalność w słabej -- zależność od ety niezbyt interesująca? Pierwsza liczba powinna być większa żeby było cokolwiek różne niż zero
    "SCF_THR": "1.e-40",
    "SCF_MAXITER": "100",
    "SCF_DIIS_START": "0",
    "SCF_DIIS_SIZE": "6",
    "CC_THR": "1.e-12",
    "CC_MAXITER": "100",
    "CC_DIIS_START": "1",
    "CC_DIIS_SIZE": "6",
    "PRINT": "1",
}


for c1 in np.arange(0.00, 1.00, 0.05):
    c1 = two_digits(c1)

    for c0 in [1.00]:
        c0 = two_digits(c0)

        for alpha in np.arange(0.00, 1.00, 0.05):
            alpha = two_digits(alpha)

            for gfac in [1.00]:
                gfac = two_digits(gfac)

                for nbas in [14]:
                    for nocc in [1]:
                        for omega in range(10, 31):
                            example_input_file["alpha"] = str(alpha)
                            example_input_file["nocc"] = str(nocc)
                            example_input_file["nbas"] = str(nbas)
                            example_input_file["omega"] = str(omega)
                            example_input_file["gfac"] = str(gfac)

                            with open("../HERMITE.INP", "w") as outfile:
                                outfile.write("\n".join(example_input_file.values()))

                            file_name = construct_file_name(nocc, nbas, omega, gfac, alpha, c0, c1)
                            printout_name = construct_file_name(nocc, nbas, omega, gfac, alpha, c0, c1, bold=True)

                            print(current_time(), printout_name + "...", end="")

                            file_path = "raw/" + file_name

                            if not os.path.isfile(file_path):
                                home = os.environ.get("HOME")
                                test_path = home + "/drake-mpfun"
                                result = subprocess.run(f"cd {test_path} && ./test", shell=True, stdout=subprocess.PIPE)

                                with open(file_path, "w") as outfile:
                                    outfile.write(result.stdout.decode())

                                output_converted.convert_output(file_name)

                            print(termcolor.colored("done", "green"))

                            if datetime.datetime.now().timestamp() % 1000 < 100:
                                time.sleep(60)
