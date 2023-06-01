# SPDX-FileCopyrightText: 2023 NORCE
# SPDX-License-Identifier: MIT

CURRENT_DIRECTORY="$PWD"

# Dune modules
for module in common geometry grid istl
do   git clone https://gitlab.dune-project.org/core/dune-$module.git --branch v2.9.0
done
for module in common geometry grid istl
do   ./dune-common/bin/dunecontrol --only=dune-$module cmake -DCMAKE_DISABLE_FIND_PACKAGE_MPI=1
     ./dune-common/bin/dunecontrol --only=dune-$module make -j5
done

# OPM modules
for repo in common grid models simulators
do
    git clone https://github.com/OPM/opm-$repo.git
done

mkdir build

for repo in common grid models
do
    mkdir build/opm-$repo
    cd build/opm-$repo
    cmake -DUSE_MPI=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CURRENT_DIRECTORY/dune-common/build-cmake;$CURRENT_DIRECTORY/dune-grid/build-cmake;$CURRENT_DIRECTORY/dune-geometry/build-cmake;$CURRENT_DIRECTORY/dune-istl/build-cmake;$CURRENT_DIRECTORY/build/opm-common;$CURRENT_DIRECTORY/build/opm-grid" $CURRENT_DIRECTORY/opm-$repo
    make -j5
    cd ../..
done    

mkdir build/opm-simulators
cd build/opm-simulators
cmake -DUSE_MPI=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CURRENT_DIRECTORY/dune-common/build-cmake;$CURRENT_DIRECTORY/dune-grid/build-cmake;$CURRENT_DIRECTORY/dune-geometry/build-cmake;$CURRENT_DIRECTORY/dune-istl/build-cmake;$CURRENT_DIRECTORY/build/opm-common;$CURRENT_DIRECTORY/build/opm-grid;$CURRENT_DIRECTORY/build/opm-models" $CURRENT_DIRECTORY/opm-simulators
make -j5 flow
cd ../..