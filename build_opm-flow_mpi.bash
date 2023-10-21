CURRENT_DIRECTORY="$PWD"

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
    cmake -DUSE_MPI=1 -DOPM_ENABLE_PYTHON=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CURRENT_DIRECTORY/build/opm-common;$CURRENT_DIRECTORY/build/opm-grid" $CURRENT_DIRECTORY/opm-$repo
    make -j5
    cd ../..
done    

mkdir build/opm-simulators
cd build/opm-simulators
cmake -DUSE_MPI=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$CURRENT_DIRECTORY/build/opm-common;$CURRENT_DIRECTORY/build/opm-grid;$CURRENT_DIRECTORY/build/opm-models" $CURRENT_DIRECTORY/opm-simulators
make -j5 flow
cd ../..