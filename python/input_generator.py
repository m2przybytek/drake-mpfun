import datetime
import os
import subprocess
import time

import output_converted

import numpy as np
import termcolor


home = os.environ.get("HOME")
test_path = home + "/drake-mpfun"


def two_digits(float_number):
    return "%.10f" % float_number


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


def update_input_file(input_file, **kwargs):
    for key, value in kwargs.items():
        input_file[key] = str(value)

    return input_file


def perform_calculation(omega=None, c1=None, c0=None, alpha=None, gfac=None, nbas=None, nocc=None, overwrite=False):
    sleep()

    c1 = two_digits(c1)
    c0 = two_digits(c0)
    alpha = two_digits(alpha)
    gfac = two_digits(gfac)

    example_input_file = {
        "calc_type": "FCCD",
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

    input_file = update_input_file(
        example_input_file, omega=omega, c1=c1, c0=c0, alpha=alpha, gfac=gfac, nbas=nbas, nocc=nocc
    )

    with open(f"{test_path}/HERMITE.INP", "w") as outfile:
        outfile.write("\n".join(input_file.values()))

    file_name = construct_file_name(nocc, nbas, omega, gfac, alpha, c0, c1)
    printout_name = construct_file_name(nocc, nbas, omega, gfac, alpha, c0, c1, bold=True)

    print(current_time(), printout_name + "...", end="")

    file_path = "raw/" + file_name

    if overwrite or not os.path.isfile(file_path):
        result = subprocess.run(f"cd {test_path} && ./test", shell=True, stdout=subprocess.PIPE)

        with open(file_path, "w") as outfile:
            outfile.write(result.stdout.decode())
    else:
        with open(file_path, "r") as infile:
            result = infile.read()

    scf, mp2, mp2ec, fccd = output_converted.convert_output(file_name)

    print(termcolor.colored("done", "green"))

    return scf, mp2, mp2ec, fccd


def sleep():
    if (datetime.datetime.now().timestamp() + 1000) % 5000 < 1000:
        time.sleep(500)


if __name__ == "__main__":
    omega_range = range(1, 2)
    c1_range = np.arange(0.00, 0.05, 1.01)
    c0_range = [1.00]
    alpha_range = np.arange(0.00, 0.05, 1.01)
    gfac_range = [1.00]
    nbas_range = [14]
    nocc_range = [1]

    for omega in omega_range:
        for c1 in c1_range:
            for c0 in c0_range:
                for alpha in alpha_range:
                    for gfac in gfac_range:
                        for nbas in nbas_range:
                            for nocc in nocc_range:
                                perform_calculation(omega, c1, c0, alpha, gfac, nbas, nocc)
