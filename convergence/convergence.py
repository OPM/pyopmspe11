# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Script to reproduce the results in https://doi.org/10.2118/231853-PA"""

import os
import csv
import subprocess
import numpy as np
from mako.template import Template


def run_command(command, wait=True):
    """Use process to execute the terminal runs"""
    with subprocess.Popen(command, shell=True) as process:
        if wait:
            process.wait()


DOMAIN = "full"  # Simulate the localized 'lower' or 'full' domain
SIZES = ["40", "20", "10", "5"]  # Full domain, example for the online docs
# SIZES = ["320", "160", "80", "40", "20", "10", "5"] # Lower domain, example for the online docs
# SIZES = ["40", "20", "10", "5", "2_5", "1_25", "62_5c"] # Simulate more refinements (full domain)
# Simulate more refinements (lower domain)
# SIZES = ["320", "160", "80", "40", "20", "10", "5", "2_5", "1_25", "62_5c", "31_25c"]
QUANTITY = (
    "tcpu"  # Summary variable for the time series from the .SMSPEC simulation files
)
COLORS = (
    "#f2d4ff,#e5a9ff,#d97eff,#cc53ff,#bf28ff,#b300ff,#9900d6,#7f00ad,#660084,#4c005b,"
    + "#330033,#1a001a"
)
NODE_FLAG = "-n lower" if DOMAIN == "lower" else ""
TIMES = "0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,350,400,450,500"
template = Template(filename="spe11b.mako")

# Run simulations
for index, size in enumerate(SIZES):
    variables = {"i": index, "domain": DOMAIN}
    rendered = template.render(**variables)
    config_name = f"{DOMAIN}_cp{index}-z{size}mish-x{size}m.toml"
    output_name = f"{DOMAIN}_cp{index}-z{size}mish-x{size}m"
    with open(config_name, "w", encoding="utf8") as file:
        file.write(rendered)
    run_command(
        f"pyopmspe11 -i {config_name} -o {output_name} {NODE_FLAG} -f 0 -w 0.1 -t {TIMES} "
        f"-r 840,1,120 -m deck_flow_data"
    )

# Generate CSV data
parallel_command = ""
for index, size in enumerate(SIZES):
    name = f"{DOMAIN}_cp{index}-z{size}mish-x{size}m"
    parallel_command += (
        f"pyopmspe11 -i {name}.toml -o {name} {NODE_FLAG} -f 0 -w 0.1 "
        + f"-t {TIMES}  -r 840,1,120 -m data -g dense_performance-spatial & "
    )
run_command(parallel_command + "wait")

# Spatial maps into a single figure (from simulation files .UNRST)
names = ""
legend = ""
for index, size in enumerate(SIZES):
    case_name = f"{DOMAIN}_cp{index}-z{size}mish-x{size}m"
    names += f"{case_name}/{case_name.upper()} "
    legend += f"{size} m  "
vals = int(np.ceil(len(SIZES) / 2))
cols = int(2.5 * vals)
run_command(
    f"plopm -i '{names.strip()}' -save {DOMAIN}_spatial_map_all -z 0 -v xco2l "
    f"-c cet_CET_CBTL1_r -clabel 'mass fraction of CO2 in liquid [-]' -d 16,{cols} "
    f"-cformat .2e -subfigs {vals},2 -cbsfax 0.35,0.97,0.3,0.02 -delax 1 -suptitle 0 "
    f"-t '{legend.strip()}'"
)

# Individual spatial maps (from simulation files .UNRST)
parallel_command = ""
for index, size in enumerate(SIZES):
    case_name = f"{DOMAIN}_cp{index}-z{size}mish-x{size}m"
    parallel_command += (
        f"plopm -i {case_name}/{case_name.upper()} "
        + f"-save {DOMAIN}_spatial_map_{case_name} -z 0 -v xco2l "
        + "-clabel 'mass fraction of CO2 in liquid [-]' -remove 1,1,1,1 "
        + "-c cet_CET_CBTL1_r -d 16,5 -cformat .2e -b '[0,6.36e-2]' & "
    )
run_command(parallel_command + "wait")

# Example of plotting a time series variable (from simulation files .SMSPEC)
run_command(
    f"plopm -i '{names.strip()}' -v {QUANTITY} -save {DOMAIN}_{QUANTITY} -ylabel "
    f"'{QUANTITY} [s]' -tunits y -labels '{legend.strip()}' -xformat .0f -e solid -lw 4"
)

# Time series from CSV files
labels = legend.split("  ")
for case_type in ["time_series", "performance_time_series_detailed"]:
    input_names = ""
    for index, size in enumerate(SIZES):
        base = f"{DOMAIN}_cp{index}-z{size}mish-x{size}m/spe11b_{case_type}"
        if index == 0:
            with open(f"{base}.csv", "r", encoding="utf8") as file:
                quantities = next(csv.reader(file))
            quantities = [q.strip() for q in quantities]
        input_names += f"{base} "
    for column_index in range(1, len(quantities)):
        case_legend = legend
        if case_type == "performance_time_series_detailed":
            case_legend = ""
            for index, size in enumerate(SIZES):
                base = f"{DOMAIN}_cp{index}-z{size}mish-x{size}m/spe11b_{case_type}"
                rows: list[list[float]] = []
                with open(f"{base}.csv", "r", encoding="utf8") as file:
                    for row_index, row in enumerate(csv.reader(file)):
                        if row_index > 1:
                            rows.append([float(value) for value in row])

                values = np.asarray(rows, dtype=float)
                if column_index in [3, 4]:
                    case_legend += (
                        labels[index] + f" max={values[:, column_index].max():.3e}  "
                    )
                else:
                    case_legend += (
                        labels[index] + f" sum={values[:, column_index].sum():.3e}  "
                    )
        variable = quantities[column_index].split(" ")[0]
        run_command(
            f"plopm -i '{input_names.strip()}' -c '{COLORS}' -csv '1,{column_index + 1}' "
            f"-save {DOMAIN}_{variable} -ylabel '{quantities[column_index]}' "
            f"-tunits y -labels '{case_legend.strip()}' -xformat .0f -e solid -lw 4"
        )

# Get the SPE11B benchmark data from the participants to compare our results
if not os.path.isdir("spe11b"):  # Skip if the data has been downloaded
    keep_files = [
        "spe11b_time_series.csv",
        "spe11b_spatial_map_500y.csv",
    ]  # To save memmory, we only keep the used files in this example
    file_ids = [
        375695,
        375747,
        375749,
        375697,
        375737,
        375744,
        375727,
        375736,
        375734,
        375753,
        375739,
        375733,
        375743,
        375740,
        375725,
        375754,
        375726,
        375745,
        375738,
        375721,
        375741,
        375751,
        375748,
        375723,
        375750,
        375722,
        375746,
        375732,
        375752,
        375742,
        375735,
        375694,
        375724,
        375698,
    ]
    for file_id in file_ids:
        run_command(
            f"curl -L -O https://darus.uni-stuttgart.de/api/access/datafile/{file_id}"
        )
        run_command(f"unzip {file_id}")
        os.remove(str(file_id))
        latest = sorted(
            name
            for name in os.listdir("spe11b")
            if os.path.isdir(os.path.join("spe11b", name))
        )[-1]
        folder_path = f"spe11b/{latest}"
        for item in os.listdir(folder_path):
            item_path = os.path.join(folder_path, item)
            if os.path.isfile(item_path) and item not in keep_files:
                os.remove(item_path)

# Example of plotting a time series variable adding all SPE11B participants (from csvs)
participant_names = ""
participant_colors = ""
for participant in sorted(
    name for name in os.listdir("spe11b") if os.path.isdir(os.path.join("spe11b", name))
):
    participant_names += f"spe11b/{participant}/spe11b_time_series "
    participant_colors += "#a8d8e3,"
for index, size in enumerate(SIZES):
    participant_names += f"{DOMAIN}_cp{index}-z{size}mish-x{size}m/spe11b_time_series "
    participant_colors += f"{COLORS.split(',')[index]},"
run_command(
    f"plopm -i '{participant_names.strip()}' -csv '1,3' -c '{participant_colors.strip(',')}' "
    f"-save {DOMAIN}_mobA_adding_participants -y '[2.1e7,4e7]' -ylabel 'mobA [kg]' -tunits y "
    "-loc empty -xlog 1 -e solid -lw 2"
)

# Spatial maps into a single figure adding all SPE11B participants (from csvs)
map_names = ""
map_titles = ""
for participant in sorted(
    name for name in os.listdir("spe11b") if os.path.isdir(os.path.join("spe11b", name))
):
    map_names += f"spe11b/{participant}/spe11b_spatial_map_500y "
    map_titles += f"{participant}  "
for index, size in enumerate(SIZES):
    map_names += f"{DOMAIN}_cp{index}-z{size}mish-x{size}m/spe11b_spatial_map_500y "
    map_titles += f"{size} m  "
grid_x = int(8 + np.ceil((len(SIZES) - 1) / 5))
grid_y = int(15.5 + 2.5 * np.floor(len(SIZES) / 2))
run_command(
    f"plopm -i '{map_names.strip()}' -csv '1,2,9' -save {DOMAIN}_spatial_map_adding_participants "
    "-z 0 -f 24 -b '[0,5e3]' -clabel 'total CO$_2$ mass [kg] at 500 years' -c cet_CET_CBTL1_r -d "
    f"60,{grid_y} -cformat .2e -subfigs {grid_x},5 -cbsfax 0.35,0.97,0.3,0.02 -delax 1 "
    f"-suptitle 0 -t '{map_titles.strip()}'"
)
