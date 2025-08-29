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
without_files = sort(glob("*without_prep*.xlsx", "./RESULTS"))
output_file = "./RESULTS/Combined_results/all_results_without_prep.xlsx"
combine_results(without_files, output_file)

with_files = sort(glob("*with_prep*.xlsx", "./RESULTS"))
output_file = "./RESULTS/Combined_results/all_results_with_prep.xlsx"
combine_results(with_files, output_file)

