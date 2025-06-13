import os
import pandas as pd
import csv

text_file_names = ["bran_sp.txt", "E10_exp2.txt", "E17_e15.txt", "E1_E6_exp2.txt", "E20.txt",
                   "E21.txt", "E22.txt", "E2.txt", "E3_exp2.txt", "E4_exp2_sp.txt",
                   "E5.txt", "est_E18.txt", "leit_spain.txt", "ord_e19_sp.txt", "pajares.txt",
                   "pena_e9.txt", "PORT.txt", "spain_pyrenes.txt", "vega_spa.txt", "yecla.txt",
                   "nor1.txt", "nor2.txt", "nor4.txt", "nor5.txt", "nor6.txt",
                   "nor7.txt", "par_nor6.txt", "s1.txt", "s2.txt", "S3.txt",
                   "S5.txt", "sve2.txt", "sve3.txt", "sve4.txt", "iceland.txt", "austria.txt",
                   "bri_fr.txt", "ch.txt", "cr_fr.txt", "D1_D2.txt", "D3.txt",
                   "D4.txt", "Dother.txt", "Fgal1.txt", "fr1_5.txt", "fr4.txt",
                   "fr6.txt", "fr_P.txt", "galiber_fr.txt", "ganon_fr.txt", "gv_fr.txt",
                   "poland.txt", "vallon_fr.txt", "italymixed.txt"]

# List of desired column names
desired_column_names = [
    "bio_1_Annual_Mean_Temperature", "bio_2_Mean_Diurnal_Range", "bio_3_Isothermality",
    "bio_4_Temperature_Seasonality", "bio_5_Max_Temperature_of_Warmest_Month",
    "bio_6_Min_Temperature_of_Coldest_Month", "bio_7_Temperature_Annual_Range",
    "bio_8_Mean_Temperature_of_Wettest_Quarter", "bio_9_Mean_Temperature_of_Driest_Quarter",
    "bio_10_Mean_Temperature_of_Warmest_Quarter", "bio_11_Mean_Temperature_of_Coldest_Quarter",
    "bio_12_Annual_Precipitation", "bio_13_Precipitation_of_Wettest_Month",
    "bio_14_Precipitation_of_Driest_Month", "bio_15_Precipitation_Seasonality",
    "bio_16_Precipitation_of_Wettest_Quarter", "bio_17_Precipitation_of_Driest_Quarter",
    "bio_18_Precipitation_of_Warmest_Quarter", "bio_19_Precipitation_of_Coldest_Quarter",
    "tavg_Jan", "tavg_Feb", "tavg_Mar", "tavg_Apr", "tavg_May", "tavg_Jun", "tavg_Jul", "tavg_Aug", "tavg_Sep", "tavg_Oct", "tavg_Nov", "tavg_Dec",
    "prec_Jan", "prec_Feb", "prec_Mar", "prec_Apr", "prec_May", "prec_Jun", "prec_Jul", "prec_Aug", "prec_Sep", "prec_Oct", "prec_Nov", "prec_Dec",
    "tmax_Jan", "tmax_Feb", "tmax_Mar", "tmax_Apr", "tmax_May", "tmax_Jun", "tmax_Jul", "tmax_Aug", "tmax_Sep", "tmax_Oct", "tmax_Nov", "tmax_Dec",
    "tmin_Jan", "tmin_Feb", "tmin_Mar", "tmin_Apr", "tmin_May", "tmin_Jun", "tmin_Jul", "tmin_Aug", "tmin_Sep", "tmin_Oct", "tmin_Nov", "tmin_Dec",
    "ai_Jan", "ai_Feb", "ai_Mar", "ai_Apr", "ai_May", "ai_Jun", "ai_Jul", "ai_Aug", "ai_Sep", "ai_Oct", "ai_Nov", "ai_Dec",
    "srad_Jan", "srad_Feb", "srad_Mar", "srad_Apr", "srad_May", "srad_Jun", "srad_Jul", "srad_Aug", "srad_Sep", "srad_Oct", "srad_Nov", "srad_Dec",
    "vapr_Jan", "vapr_Feb", "vapr_Mar", "vapr_Apr", "vapr_May", "vapr_Jun", "vapr_Jul", "vapr_Aug", "vapr_Sep", "vapr_Oct", "vapr_Nov", "vapr_Dec",
    "wind_Jan", "wind_Feb", "wind_Mar", "wind_Apr", "wind_May", "wind_Jun", "wind_Jul", "wind_Aug", "wind_Sep", "wind_Oct", "wind_Nov", "wind_Dec", "Elev"]

excel_file_path = r"C:\Users\sayeraselvan\Desktop\climate_data.xlsx"

output_directory = r"C:\Users\sayeraselvan\Desktop\txt1"
os.makedirs(output_directory, exist_ok=True)

excel_data = pd.read_excel(excel_file_path)

for desired_column_name in desired_column_names:
    results = {}

    for text_file_name in text_file_names:
        text_file_path = os.path.join(r"C:\Users\sayeraselvan\Desktop\txt\all", text_file_name)
        with open(text_file_path, "r") as file:
            text_values = [int(line.strip()) if line.strip().isdigit() else None for line in file]

        corresponding_cell_values = []
        non_matching_text_values = []

        for text_value in text_values:
            matching_index = excel_data[excel_data.iloc[:, 1] == text_value].index

            if not matching_index.empty:
                corres_value = excel_data.at[matching_index[0], desired_column_name]
                corresponding_cell_values.append(corres_value)
            else:
                corresponding_cell_values.append(None)
                non_matching_text_values.append(text_value)

        filtered_values = [value for value in corresponding_cell_values if value is not None]

        sum_of_values = sum(filtered_values)
        mean_value = round(sum_of_values / len(filtered_values), 2)

        results[text_file_name] = mean_value

    # Create a filename for the CSV file
    output_csv_file = os.path.join(output_directory, f"results_{desired_column_name}.csv")

    csv_data = [{"Filename": filename.split(".")[0] + "v", desired_column_name: ai_apr} for filename, ai_apr in results.items()]

    with open(output_csv_file, "w", newline="") as csv_file:
        fieldnames = ["Filename", desired_column_name]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)

        writer.writeheader()  # Write the header row

        for row in csv_data:
            writer.writerow(row)

    print(f"Data for '{desired_column_name}' has been saved to {output_csv_file}")
