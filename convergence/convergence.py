# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: MIT

import os
import csv
import numpy as np
from mako.template import Template

domain = "full" # Simulate the localized 'lower' or 'full' domain
sizes = ["40", "20", "10", "5"] # Full domain, example for the online docs
# sizes = ["320", "160", "80", "40", "20", "10", "5"] # Lower domain, example for the online docs
# sizes = ["40", "20", "10", "5", "2_5", "1_25", "62_5c"] # Simulate more refinements (full domain)
# sizes = ["320", "160", "80", "40", "20", "10", "5", "2_5", "1_25", "62_5c", "31_25c"] # Simulate more refinements (lower domain)

# Run the simulations
nflag = "-n lower" if domain == "lower" else ""
mytemplate = Template(filename="spe11b.mako")
for i, size in enumerate(sizes):
    var = {"i": i, "domain": domain}
    filledtemplate = mytemplate.render(**var)
    with open(
        f"{domain}_cp{i}-z{size}mish-x{size}m.toml",
        "w",
        encoding="utf8",
    ) as file:
        file.write(filledtemplate)
    os.system(f"pyopmspe11 -i {domain}_cp{i}-z{size}mish-x{size}m.toml -o {domain}_cp{i}-z{size}mish-x{size}m {nflag} -f 0 -w 0.1 -t 0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,350,400,450,500 -r 840,1,120 -m deck_flow_data")

# Generate the csv data
command = ""
for i, size in enumerate(sizes):
    command += f"pyopmspe11 -i {domain}_cp{i}-z{size}mish-x{size}m.toml -o {domain}_cp{i}-z{size}mish-x{size}m {nflag} -f 0 -w 0.1 -t 0,5,10,15,20,25,30,35,40,45,50,75,100,150,200,250,300,350,400,450,500 -r 840,1,120 -m data -g dense_performance-spatial & "
command += "wait"
os.system(command)

# Spatial maps into a single figure (from simulation files .UNRST)
names, legend = "", ""
for i, size in enumerate(sizes):
    name = f"{domain}_cp{i}-z{size}mish-x{size}m"
    names += f"{name}/{name.upper()} "
    legend += f"{size} m  "
os.system(f"plopm -i '{names[:-1]}' -save {domain}_spatial_map_all -z 0 -v xco2l -c cet_CET_CBTL1_r -clabel 'mass fraction of CO2 in liquid [-]' -d 16,{int(2.5*np.ceil(len(sizes)/2))} -cformat .2e -subfigs {int(np.ceil(len(sizes)/2))},2 -cbsfax 0.35,0.97,0.3,0.02 -delax 1 -suptitle 0 -t '{legend}'")

# Individual spatial maps (from simulation files .UNRST)
command = ""
for i, size in enumerate(sizes):
    name = f"{domain}_cp{i}-z{size}mish-x{size}m"
    command += f"plopm -i {name}/{name.upper()} -save {domain}_spatial_map_{name} -z 0 -v xco2l -clabel 'mass fraction of CO2 in liquid [-]' -remove 1,1,1,1 -c cet_CET_CBTL1_r -d 16,5 -cformat .2e -b '[0,6.36e-2]' & "
command += "wait"
os.system(command)

# Example of plotting a time series variable (from simulation files .SMSPEC)
colors = '#f2d4ff,#e5a9ff,#d97eff,#cc53ff,#bf28ff,#b300ff,#9900d6,#7f00ad,#660084,#4c005b,#330033,#1a001a'
quan = "tcpu"
os.system(f"plopm -i '{names[:-1]}' -v {quan} -save {domain}_{quan} -ylabel '{quan} [s]' -tunits y -labels '{legend[:-2]}' -xformat .0f -e solid -lw 4")

# Time series (from csvs)
legends = legend.split("  ")
for case in ["time_series", "performance_time_series_detailed"]:
    names = ""
    for i, size in enumerate(sizes):
        name = f"{domain}_cp{i}-z{size}mish-x{size}m/spe11b_{case}"
        if i == 0:
            with open(
                f"{name}.csv",
                "r",
                encoding="utf8",
            ) as file:
                for j, row in enumerate(csv.reader(file)):
                    quans = row
                    break
            quans = [quan.strip() for quan in quans]
        names += f"{name} "

    for n in range(1,len(quans)):
        if case == "performance_time_series_detailed":
            legend = ""
            for i, size in enumerate(sizes):
                name = f"{domain}_cp{i}-z{size}mish-x{size}m/spe11b_{case}"
                csvs = []
                with open(
                    f"{name}.csv",
                    "r",
                    encoding="utf8",
                ) as file:
                    for j, row in enumerate(csv.reader(file)):
                        if j > 1:
                            csvs.append(
                                [float(column) for column in row]
                            )
                csvs = np.array(csvs)
                if n in [3, 4]:
                    legend += legends[i] + f" max={max(csvs[:,n]):.3e}  "
                else:
                    legend += legends[i] + f" sum={sum(csvs[:,n]):.3e}  "
        var = quans[n].split(" ")[0]
        os.system(f"plopm -i '{names[:-1]}' -c '{colors}' -csv '1,{n+1}' -save {domain}_{var} -ylabel '{quans[n]}' -tunits y -labels '{legend[:-2]}' -xformat .0f -e solid -lw 4")

# Get the SPE11B benchmark data from the participants to compare our results
if not os.path.isdir("spe11b"): # Skip if the data has been downloaded
    keep = ["spe11b_time_series.csv", "spe11b_spatial_map_500y.csv"] # To save memmory, we only keep the used files in this example
    for id in [375695,375747,375749,375697,375737,375744,375727,375736,375734,375753,375739,375733,375743,375740,375725,375754,375726,375745,375738,375721,375741,375751,375748,375723,375750,375722,375746,375732,375752,375742,375735,375694,375724,375698]:
        os.system(f"curl -L -O https://darus.uni-stuttgart.de/api/access/datafile/{id}")
        os.system(f"unzip {id}")
        os.remove(str(id))
        folder = sorted(name for name in os.listdir("spe11b") if os.path.isdir(os.path.join("spe11b", name)))[-1]
        folderp = f"spe11b/{folder}"
        for item in os.listdir(folderp):
            itemp = os.path.join(folderp, item)
            if os.path.isfile(itemp) and item not in keep:
                os.remove(itemp)

# Example of plotting a time series variable adding all SPE11B participants (from csvs)
names, colours = "", ""
for part in sorted(name for name in os.listdir("spe11b") if os.path.isdir(os.path.join("spe11b", name))):
    names += f"spe11b/{part}/spe11b_time_series "
    colours += '#a8d8e3,'
for i, size in enumerate(sizes):
    names += f"{domain}_cp{i}-z{size}mish-x{size}m/spe11b_time_series "
    colours += f"{colors.split(',')[i]},"
os.system(f"plopm -i '{names[:-1]}' -csv '1,3' -c '{colours[:-1]}' -save {domain}_mobA_adding_participants -y '[2.1e7,4e7]' -ylabel 'mobA [kg]' -tunits y -loc empty -xlog 1 -e solid -lw 2")

# Spatial maps into a single figure adding all SPE11B participants (from csvs)
names, titles = "", ""
for part in sorted(name for name in os.listdir("spe11b") if os.path.isdir(os.path.join("spe11b", name))):
    names += f"spe11b/{part}/spe11b_spatial_map_500y "
    titles += f"{part}  "
for i, size in enumerate(sizes):
    names += f"{domain}_cp{i}-z{size}mish-x{size}m/spe11b_spatial_map_500y "
    titles += f"{size} m  "
os.system(f"plopm -i '{names[:-1]}' -csv '1,2,9' -save {domain}_spatial_map_adding_participants -z 0 -f 24 -b '[0,5e3]' -clabel 'total CO$_2$ mass [kg] at 500 years' -c cet_CET_CBTL1_r -d 60,{15.5 + 2.5*np.floor(len(sizes)/2)} -cformat .2e -subfigs {int(8+np.ceil((len(sizes)-1)/5))},5 -cbsfax 0.35,0.97,0.3,0.02 -delax 1 -suptitle 0 -t '{titles[:-1]}'")
