. tests/scripts/initialize_output_folders.sh &
. tests/scripts/get_plopm.sh &
wait

. tests/scripts/docs_hello_world.sh &
. tests/scripts/docs_cp_grids.sh &
. tests/scripts/docs_localized_lower_domain.sh &
wait

. tests/scripts/docs_check_outputs.sh
