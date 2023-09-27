
# Ionosonde-Analysis-Toolkit

Ionosonde Analysis Toolkit (IAT), consisting of four programs to make it easier for researchers to interact effectively with non-digisonde systems like CADI  

* mdx2txt_cli.py

This Python program, named MDX2TXT Converter, is designed to convert MD2 and MD4 ionogram data files into human-readable TXT format. It provides a command-line interface for easy batch processing of MD2 and MD4 files in a specified input directory and saves the converted files in an output directory.  
Example Usage:

To convert a single MD2 file located at input_file.md2 and save the result in the directory output_folder, run:

    python mdx2txt.py input_file.md2 output_folder

To convert all MD2 and MD4 files in a folder called data_folder and save the results in a directory named output_data, run:

    python mdx2txt.py data_folder output_data
