import os
import subprocess
import time

from ase.parallel import parprint, paropen

from hl_gaps_pub.hl_gaps_pub import (
    calculate_gap,
    embed_confs,
    parse_sdf_db,
)


def main():
    """
    Main function to process molecules from a database, calculate their HOMO-LUMO gaps,
    and save the results to a file.
    """

    # Environment variable for xTB
    os.environ["XTBPATH"] = "/home/sat/miniforge3/envs/py310hl_gaps_pub/bin/xtb"

    # Input parameters
    nconfs = 1  # Number of conformers to generate
    accuracy = 1.0  # Accuracy of calculations (smaller means higher accuracy, e.g., 0.0001)
    temp = 300.0  # Electronic temperature
    method = "GFN2-xTB"  # Tight-binding method to use ["GFN0-xTB","GFN1-xTB","GFN2-xTB","IPEA-xTB"]
    id_start = 407268
    id_end = 407270

    # Initialize results file
    with paropen("results.raw", mode="a") as ff:
        ff.write("#   ID   GAP    TIME SMILE \n")

    t0 = time.time()

    for ii in range(id_start, id_end):
        t0i = time.time()
        pref = subprocess.getoutput("pwd") + "/../data/coconut/db_split"
        try:
            coconut_df = parse_sdf_db(f"{pref}/{ii}_COCONUT_2022_01_2D.SDF")
        except Exception as e:
            print(f"Error parsing SDF for ID {ii}: {e}")
            continue
        smilestr = coconut_df["SMILES"][0]
        parprint(ii, smilestr)
        confs = embed_confs(smile=smilestr, num_confs=nconfs)
        gap = calculate_gap(
            molecule=confs, method=method, accuracy=accuracy, temperature=temp
        )

        t1i = time.time()
        res_l = [ii, gap, t1i - t0i, smilestr]
        parprint("finished molecule: ", res_l)
        with paropen("results.raw", mode="a") as ff:
            ff.write(
                "{:>6d} {:>5.2f} {:>7.1f} {:s}\n".format(
                    res_l[0], res_l[1], res_l[2], res_l[3]
                )
            )

    t1 = time.time()
    print(f"Total runtime {t1 -t0} seconds")

if __name__ == "__main__":
    main()
