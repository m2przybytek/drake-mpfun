import os

import input_generator


def optimize(param1, param2, **other_params):
    cutoff = other_params.pop("cutoff")

    param2_key = list(param2)[0]
    param2_value = param2[param2_key][0]

    for param, param_specs in param1.items():
        start, delta = param_specs
        end = start + delta
        if end < 0:
            start, end = end, start
            delta = -delta

        _, _, _, energy_start = input_generator.perform_calculation(
            **{param: start, param2_key: param2_value, **other_params}
        )
        print(param, start, param2_key, param2_value, energy_start)

        _, _, _, energy_end = input_generator.perform_calculation(
            **{param: end, param2_key: param2_value, **other_params}
        )
        print(param, end, param2_key, param2_value, energy_end)

        if abs(energy_end - energy_start) < cutoff:
            return (energy_start + energy_end) / 2, param, end, param2_key, param2_value

        if energy_end < energy_start:
            return optimize({param: (end, delta)}, param2, cutoff=cutoff, **other_params)

        delta = -delta / 2

        return optimize(param2, {param: (end, delta)}, cutoff=cutoff, **other_params)


if __name__ == "__main__":
    for omega in range(25, 35, 1):
        cutoff = 0.000000000001
        c0 = 1.00
        gfac = 1.00
        nbas = 14
        nocc = 1

        other_params = {"omega": omega, "c0": c0, "gfac": gfac, "nbas": nbas, "nocc": nocc}

        result, param1, param1_value, param2, param2_value = optimize(
            {"c1": [0.20, 0.05]}, {"alpha": [0.00, 0.05]}, cutoff=cutoff, **other_params
        )
        print(
            f"Optimized {param1} to be {param1_value} and {param2} to be {param2_value} with energy {result}"
            f"for other params equal {other_params}"
        )

        omega = input_generator.two_digits(omega)
        c0 = input_generator.two_digits(c0)
        gfac = input_generator.two_digits(gfac)
        nbas = input_generator.two_digits(nbas)
        nocc = input_generator.two_digits(nocc)

        with open(f"optimization/omega={omega}c0={c0}gfac={gfac}nbas={nbas}cutoff={cutoff}", "w") as outfile:
            outfile.write(str(result))
            outfile.write("\n")
            outfile.write(str(param1))
            outfile.write(" ")
            outfile.write(str(param1_value))
            outfile.write("\n")
            outfile.write(str(param2))
            outfile.write(" ")
            outfile.write(str(param2_value))
