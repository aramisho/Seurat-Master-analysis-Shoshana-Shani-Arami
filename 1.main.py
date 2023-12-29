import pandas as pd
import fileReader
from fileReader import FilesInPath

#incorporated and matched relevant cell information from the metadata file, including labels such as class, cluster, subclass, and region

metadata = pd.read_csv('/home/local/BGU-USERS/aramisho/py_projects/allen_brain_mouse/metadata (1).csv')
transpose_matrix_files = FilesInPath('/home/local/BGU-USERS/aramisho/py_projects/allen_brain_mouse/matrix_split/splits')
output_path = '/home/local/BGU-USERS/aramisho/py_projects/allen_brain_mouse/matrix_split/non_normalized_by_region/'
metadata = metadata[["sample_name", "class_label", "cluster_label", "subclass_label", "region_label"]]
print('start')
# START - Split the matrix to file by Cluster_Label


for file, path in transpose_matrix_files:
    #  Merge the files.
    try:
        file.rename(columns={'Unnamed: 0': 'sample_name'}, inplace=True)
        file['sample_name'] = file['sample_name'].apply(lambda x: f'{x}'.replace('.', '-'))
        try:
            merged_table = metadata.merge(file, on='sample_name', how='inner')
        except:
            print(f'error in {path} somthing wrong with the inner join.')
    except:
        print(f'error in {path} sample_name column doesnt exists.')

    regions = merged_table['region_label'].unique()
    print(file)
    print(regions)
    # For each region in the current file join the required parameters.
    for region in regions:
        region = str(region)
        region = region.replace('/', '\\')
        reads_to_add = merged_table[merged_table['region_label'] == region]

        try:
            region_csv = pd.read_csv(output_path + region + '.csv')
            region_csv.append(reads_to_add)
        except:
            region_csv = reads_to_add

        region_csv.to_csv(output_path + region + '.csv', index=False)
# END - Split the matrix to file by Cluster_Label
# column_to_delete = ["exp_component_name", "platform_label", "cluster_color", "cluster_order", "class_color",
#                     "class_order", "subclass_color", "subclass_order", "full_genotype_color", "full_genotype_id",
#                     "full_genotype_label", "sex_color", "sex_id", "donor_sex_label", "region_color", "region_id",
#                     "cell_type_accession_color", "cell_type_accession_id", "cell_type_accession_label",
#                     "cell_type_alias_color", "cell_type_alias_id", "cell_type_alias_label", "cell_type_alt_alias_color",
#                     "cell_type_alt_alias_id", "cell_type_alt_alias_label", "cell_type_designation_color",
#                     "cell_type_designation_id", "cell_type_designation_label", "neighborhood_label", "neighborhood_id",
#                     "neighborhood_color", "external_donor_name_color", "external_donor_name_id",
#                     "external_donor_name_label", "facs_population_plan_color", "facs_population_plan_id",
#                     "facs_population_plan_label", "injection_materials_color", "injection_materials_id",
#                     "injection_materials_label", "injection_method_color", "injection_method_id",
#                     "injection_method_label", "injection_roi_color", "injection_roi_id", "injection_roi_label",
#                     "injection_type_color", "injection_type_id", "injection_type_label", "cortical_layer_label",
#                     "outlier_call", "outlier_type"]
# long_files = FilesInPath(base_region_path)
# # START - Creating Short files for the matrix files who splited by Cluster Label
# for file, path in long_files:
#     output_path = '/home/local/BGU-USERS/tzahben/shoshana/cluster_label_fixed/' + path.split('/')[1]
#     print(output_path)
#     for column in column_to_delete:
#         del file[column]
#     file.to_csv(output_path, index=False)
# print('finish')
