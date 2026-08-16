files=(
    "test_outputs/docs_localized_lower_domain/lower_domain_pvtnum-1-satnum_i,14,k_t5.png"
    "test_outputs/docs_via_deck_hello_world/spe11b_performance.png"
    "test_outputs/docs_via_deck_hello_world/spe11b_tco2_2Dmaps.png"
    "test_outputs/docs_via_deck_hello_world/spe11b_performance_detailed.png"
    "test_outputs/docs_via_deck_hello_world/spe11b_time_series_csv.png"
    "test_outputs/docs_via_deck_hello_world/spe11b_sparse_data.png"
    "test_outputs/docs_via_deck_hello_world/isothermal_sgas_i,1,k_t5.png"
    "test_outputs/docs_cp_grids/11_levels_dz_i,1,k_t5.png"
)

missing_file="test_outputs/missing_docs_files.txt"
missing=0

rm -f "$missing_file"

for f in "${files[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "$f" >> "$missing_file"
        ((missing++))
    fi
done

if (( missing == 0 )); then
    echo "All figures and files exist."
else
    echo "$missing figure(s) or file(s) missing."
    echo "See $missing_file"
fi
