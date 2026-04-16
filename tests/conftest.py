# SPDX-FileCopyrightText: 2026 NORCE Research AS
# SPDX-License-Identifier: MIT

"""Common folders for the tests"""

import pytest


@pytest.fixture(scope="session")
def test_c_i_workdir(tmp_path_factory):
    """Folder for the single spe11a run and data format for all cases"""
    return tmp_path_factory.mktemp("test_c_i_spe11a_format-spe11abc")


@pytest.fixture(scope="session")
def test_d_f_g_workdir(tmp_path_factory):
    """Folder for the single spe11b run, data, and figures generation"""
    return tmp_path_factory.mktemp("test_d_f_g_spe11b_data_plot")


@pytest.fixture(scope="session")
def test_e_h_workdir(tmp_path_factory):
    """Folder for the single spe11c run, data, and compare functionality"""
    return tmp_path_factory.mktemp("test_e_h_spe11c_compare")
