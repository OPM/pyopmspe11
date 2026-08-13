WHR="examples/spe11b"
OUT="test_outputs/docs_cp_grids"
. tests/scripts/initialize_output_folders.sh $OUT
. tests/scripts/get_plopm.sh
pyopmspe11 -i $WHR.toml -o $OUT/18_levels -f 0
cp $WHR.toml $OUT/spe11b_11-levels.toml
sed -i.bak "s/2, 2, 1, 1, 2, 1, 3, 2, 6, 2, 3, 2, 2, 8, 4, 8, 8, 1/2, 2, 2, 3, 2, 2, 8, 4, 8, 8, 1/g" $OUT/spe11b_11-levels.toml && rm -f $OUT/spe11b_11-levels.toml.bak
pyopmspe11 -i $OUT/spe11b_11-levels.toml -o $OUT/11_levels -f 0
plopm -i "$OUT/18_levels/18_LEVELS $OUT/11_levels/11_LEVELS" -o $OUT -v dz -subfigs 2,1 -delax 1 -z 0 -suptitle 0 -grid 'black,1e-2' -cbsfax 0.35,0.97,0.3,0.02
