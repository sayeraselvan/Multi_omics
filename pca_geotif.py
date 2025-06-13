import os
import rasterio
from rasterio.windows import Window

input_directory = r"C:\Users\sayeraselvan\Desktop\Visual studio"
output_directory = r"C:\Users\sayeraselvan\Desktop\Visual studio\crop"
os.makedirs(output_directory, exist_ok=True)
crop_coordinates = (25, 80, -27, 55)

# List the GeoTIFF files in the input directory.
geotiff_files = [os.path.join(input_directory, file) for file in os.listdir(input_directory) if file.endswith('.tif')]

# Function to crop and save a GeoTIFF file.
def crop_and_save_geotiff(input_tif, output_tif, crop_coordinates):
    with rasterio.open(input_tif) as src:
        left, right, bottom, top = crop_coordinates

        # Calculate the window coordinates based on the latitude and longitude coordinates.
        transform = src.transform
        col_start, col_stop = ~transform * (left, right)
        row_start, row_stop = ~transform * (top, bottom)
        window = Window(col_start, row_start, col_stop - col_start, row_stop - row_start)

        # Read and crop the data
        data = src.read(window=window)

        # Update the metadata for the output GeoTIFF
        profile = src.profile
        profile.update({
            'height': data.shape[1],
            'width': data.shape[2],
            'transform': src.window_transform(window)
        })

        # Create and save the cropped GeoTIFF
        with rasterio.open(output_tif, 'w', **profile) as dst:
            dst.write(data)

# Loop over each GeoTIFF file, crop it, and save the cropped version to the output directory.
for geotiff_file in geotiff_files:
    input_tif = geotiff_file
    output_tif = os.path.join(output_directory, os.path.basename(geotiff_file))
    crop_and_save_geotiff(input_tif, output_tif, crop_coordinates)

print("Cropping completed.")