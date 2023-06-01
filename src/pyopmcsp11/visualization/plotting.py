# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

""""
Script to plot the results
"""

from datetime import timedelta
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ecl.summary import EclSum
from ecl.eclfile import EclFile
from ecl.grid import EclGrid


def main():
    """Postprocessing"""
    parser = argparse.ArgumentParser(description="Main script to plot the results")
    parser.add_argument(
        "-t",
        "--time",
        default=0.0,
        help="The workflow time.",
    )
    parser.add_argument(
        "-p",
        "--path",
        default="output",
        help="The path to the opm simulations.",
    )
    parser.add_argument(
        "-d",
        "--deck",
        default="CSP11A",
        help="The simulated case.",
    )
    cmdargs = vars(parser.parse_known_args()[0])
    dic = {"path": cmdargs["path"].strip()}
    dic["case"] = cmdargs["deck"].strip()
    dic["rhog_ref"] = 1.86843  # CO2 reference density
    # The current implementation only computes the mass for the immiscible model
    dic["quantity"] = ["mass", "saturation", "pressure"]
    case = dic["path"] + "/output/" + f"{dic['case'].upper()}"
    dic["rst"] = EclFile(case + ".UNRST")
    dic["ini"] = EclFile(case + ".INIT")
    dic["grid"] = EclGrid(case + ".EGRID")
    dic["smsp"] = EclSum(case + ".SMSPEC")
    dic["saturation"] = np.array(dic["rst"].iget_kw("SGAS"))
    dic["pressure"] = np.array(dic["rst"].iget_kw("PRESSURE"))
    dic["density"] = np.array(dic["rst"].iget_kw("GAS_DEN"))
    dic["actnum"] = np.array(dic["grid"].export_actnum())
    dic["phiv"] = np.array(dic["ini"].iget_kw("PORV")[0])[dic["actnum"] == 1]
    if dic["case"] == "CSP11C":
        dic["dims"] = dic["grid"].get_dims()
        dic["mid_slice"] = []
        for k in range(dic["dims"][2]):
            for i in range(dic["dims"][0]):
                dic["mid_slice"].append(
                    dic["grid"].global_index(
                        ijk=(i, round(dic["grid"].get_ny() / 2), k)
                    )
                )
    for quantity in dic["quantity"]:
        dic[f"{quantity}_array"] = []
        mesh = np.array([0.0 for _ in range(dic["grid"].get_global_size())])
        for i in range(dic["rst"].num_report_steps()):
            if quantity == "mass":
                mesh[dic["actnum"] == 1] = np.array(
                    dic["saturation"][i] * dic["density"][i] * dic["phiv"]
                )
            else:
                mesh[dic["actnum"] == 1] = np.array(dic[f"{quantity}"][i])
            if dic["case"] == "CSP11C":
                dic[f"{quantity}_array"].append(
                    np.array([mesh[j] for j in dic["mid_slice"]])
                )
            else:
                dic[f"{quantity}_array"].append(mesh)
    final_time_maps(dic)
    print(f"Total simulation time: {timedelta(seconds=float(cmdargs['time']))}")


def final_time_maps(dic):
    """
    Function to plot the 2D maps for the different reservoirs and quantities

    Args:
        dic (dict): Global dictionary with required parameters

    """
    dic["boxi"] = dic["grid"].getNodePos(0, 0, 0)
    dic["boxf"] = dic["grid"].getNodePos(
        dic["grid"].getNX(),
        dic["grid"].getNY(),
        dic["grid"].getNZ(),
    )
    dic["xmx"] = np.linspace(dic["boxi"][0], dic["boxf"][0], dic["grid"].get_nx() + 1)
    dic["zmz"] = np.linspace(dic["boxi"][2], dic["boxf"][2], dic["grid"].get_nz() + 1)
    dic["xcor"], dic["zcor"] = np.meshgrid(dic["xmx"], dic["zmz"][::-1])
    for quantity in dic["quantity"]:
        dic[f"{quantity}_plot"] = np.zeros([len(dic["zmz"]) - 1, len(dic["xmx"]) - 1])
        if quantity == "mass":
            total_mass = max(dic["mass_array"][-1].sum(), 1e-12)
            for i in np.arange(0, len(dic["zmz"]) - 1):
                dic[f"{quantity}_plot"][i, :] = (
                    dic[f"{quantity}_array"][-1][
                        i * (len(dic["xmx"]) - 1) : (i + 1) * (len(dic["xmx"]) - 1)
                    ]
                    / total_mass
                )
        else:
            for i in np.arange(0, len(dic["zmz"]) - 1):
                dic[f"{quantity}_plot"][i, :] = dic[f"{quantity}_array"][-1][
                    i * (len(dic["xmx"]) - 1) : (i + 1) * (len(dic["xmx"]) - 1)
                ]
        fig, axis = plt.subplots()
        if quantity == "mass":
            imag = axis.pcolormesh(
                dic["xcor"],
                dic["zcor"],
                dic[f"{quantity}_plot"],
                shading="flat",
                cmap="jet",
            )
            axis.set_title(
                f"normalized {quantity} (Total mass: "
                + f"{dic['mass_array'][-1].sum() : .2E} [kg])"
            )
        elif quantity == "pressure":
            imag = axis.pcolormesh(
                dic["xcor"],
                dic["zcor"],
                dic[f"{quantity}_plot"],
                shading="flat",
                cmap="jet",
            )
            axis.set_title(f"{quantity} [bar]")
        else:
            imag = axis.pcolormesh(
                dic["xcor"],
                dic["zcor"],
                dic[f"{quantity}_plot"],
                shading="flat",
                cmap="jet",
            )
            axis.set_title(f"{quantity}")
        axis.axis("scaled")
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        vect = np.linspace(
            dic[f"{quantity}_plot"].min(),
            dic[f"{quantity}_plot"].max(),
            5,
            endpoint=True,
        )
        fig.colorbar(
            imag,
            cax=cax,
            orientation="vertical",
            ticks=vect,
            format=lambda x, _: f"{x:.2f}",
        )
        imag.set_clim(
            dic[f"{quantity}_plot"].min(),
            dic[f"{quantity}_plot"].max(),
        )
        fig.savefig(f"{dic['case'].upper()}_{quantity}.png", bbox_inches="tight")
        plt.close()

    return dic


if __name__ == "__main__":
    main()
