NCPUS=${1:-16}
OUT="test_outputs/paper_convergence"
. tests/scripts/initialize_output_folders.sh $OUT
. tests/scripts/get_plopm.sh
cp -r convergence/. $OUT
sed -i.bak "s/mpirun -np 8 flow --/mpirun -np $NCPUS flow --partition-method=metis --edge-weights-method=logtrans --/g" $OUT/spe11b.mako && rm -f $OUT/spe11b.mako.bak
cd $OUT
python3 convergence.py
cd ../..

files="
test_outputs/paper_convergence/full_p1.png
test_outputs/paper_convergence/full_tlinsol.png
test_outputs/paper_convergence/full_immB.png
test_outputs/paper_convergence/full_p2.png
test_outputs/paper_convergence/full_immA.png
test_outputs/paper_convergence/full_spatial_map_full_cp0-z40mish-x40m.png
test_outputs/paper_convergence/full_mobA_adding_participants.png
test_outputs/paper_convergence/full_spatial_map_full_cp2-z10mish-x10m.png
test_outputs/paper_convergence/full_spatial_map_adding_participants.png
test_outputs/paper_convergence/full_sealTot.png
test_outputs/paper_convergence/full_nliter.png
test_outputs/paper_convergence/full_nres.png
test_outputs/paper_convergence/full_liniter.png
test_outputs/paper_convergence/full_sealA.png
test_outputs/paper_convergence/full_sealB.png
test_outputs/paper_convergence/full_runtime.png
test_outputs/paper_convergence/full_spatial_map_full_cp3-z5mish-x5m.png
test_outputs/paper_convergence/full_dof.png
test_outputs/paper_convergence/full_boundTot.png
test_outputs/paper_convergence/full_mobA.png
test_outputs/paper_convergence/full_fsteps.png
test_outputs/paper_convergence/full_tcpu.png
test_outputs/paper_convergence/full_mobB.png
test_outputs/paper_convergence/full_spatial_map_full_cp1-z20mish-x20m.png
test_outputs/paper_convergence/full_MC.png
test_outputs/paper_convergence/full_tstep.png
test_outputs/paper_convergence/full_spatial_map_all.png
test_outputs/paper_convergence/full_dissB.png
test_outputs/paper_convergence/full_mass.png
test_outputs/paper_convergence/full_dissA.png
"

missing_file="test_outputs/missing_publication_files.txt"
missing=0

rm -f "$missing_file"

printf '%s\n' "$files" | while IFS= read -r f; do
    [ -z "$f" ] && continue
    if [ ! -f "$f" ]; then
        echo "$f" >> "$missing_file"
        missing=$((missing + 1))
    fi
done

if [ "$missing" -eq 0 ]; then
    echo "All figures and files exist."
    return 0
else
    echo "$missing figure(s) or file(s) missing."
    echo "See $missing_file"
    return 1
fi
