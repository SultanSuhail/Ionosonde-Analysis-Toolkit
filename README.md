# Ionosonde Analysis Toolkit (IAT)

The Ionosonde Analysis Toolkit (IAT) comprises four programs designed to simplify interactions with non-digisonde systems like CADI for researchers.

## mdx2txt_cli.py

**Description**: This program offers a command-line interface for effortless batch processing of MD2 and MD4 files within a specified input directory, saving the converted ASCII files to an output directory.

**Example Usage**:

To convert a single MD2 file located at `input_file.md2` and save the result in the directory `output_folder`, run:

```bash
python mdx2txt.py input_file.md2 output_folder
```

To convert all MD2 and MD4 files in a folder named `data_folder` and save the results in a directory named `output_data`, run:

```bash
python mdx2txt.py data_folder output_data
```

## mdx2csv_cli.py

**Description**: This Python program, named `mdx2csv_cli.py`, is a command-line utility for extracting essential features from MD2 and MD4 Ionogram data files into CSV format. The script provides flexibility in selecting which features to save in the output CSV files.

**Usage**:

To use the script, execute it from the command line with the following syntax:

```bash
python mdx2csv_cli.py [command args]
```

**Command Arguments**:

- `input_path` (str): Specifies the input folder or file to convert, which can contain one or more MD2 or MD4 ionogram data files.

- `output_directory` (str): Specifies the output directory where the converted CSV files will be saved.

- `--custom` (str, optional): Allows you to specify custom features to save using a shorthand notation. You can use the following characters to indicate which features to save:
    - `d`: Date information
    - `t`: Time information
    - `f`: Frequency information
    - `h`: Height information
    - `m`: Mean power information

**Example Usage**:

Here's an example of how to use the script:

```bash
python mdx2csv_cli.py input_folder/ output_folder/ --custom dtfhm
```

