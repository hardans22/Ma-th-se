using XLSX, DataFrames, Glob

function combine_results(file_list, output_file)
    all_data = DataFrame[]
    for file in file_list
        println(file)
        df = DataFrame(XLSX.readtable(file, "Sheet1"))
        push!(all_data, df)
    end
    combined = vcat(all_data..., cols=:union)
        
    # Sauvegarder
    XLSX.writetable(output_file, combined, overwrite = true)
end
init_path = "/home/harcenage/Téléchargements/RESULTS/"
for first_fold in ["M01_W1/", "M01_W2/", "M01_W3/", "M01_W4/"]
    for second_fold in ["MS_5/", "MS_10/"]
        without_files = sort(glob("*_without_*.xlsx", init_path*first_fold*second_fold))
        combine_results(without_files, init_path*first_fold*second_fold*"COMBINED/all_results_without_prep.xlsx")
        with_files = sort(glob("*_with_*.xlsx", init_path*first_fold*second_fold))
        combine_results(with_files, init_path*first_fold*second_fold*"COMBINED/all_results_with_prep.xlsx")
    end
end


