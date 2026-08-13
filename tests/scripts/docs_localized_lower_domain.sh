WHR="examples/spe11c"
OUT="test_outputs/docs_localized_lower_domain"
. tests/scripts/initialize_output_folders.sh $OUT
. tests/scripts/get_plopm.sh
pyopmspe11 -i $WHR.toml -o $OUT/lower_domain -f 0 -n lower
plopm -i $OUT/lower_domain/LOWER_DOMAIN -o $OUT -s ,14, -y '[1200,700]' -z 0 -grid 'black,1e-2' -t "SPE11C Cartesian lower domain (y = 2500 m)" -clabel "Facies" -c '161;163;160 101;64;147 81;124;66 181;73;57 193;127;97 127;148;191 193;147;56' -cticks '[7, 6, 5, 4, 3, 2, 1]' -v 'pvtnum - 1 - satnum'
