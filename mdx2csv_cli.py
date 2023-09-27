import argparse
import csv
import math
import os
import struct
import datetime
import numpy as np
from tqdm import tqdm
from statistics import median


# Function to convert MD2/MD4 radar data files to CSV format, extracting key features
def mdx2csv(input_file, output_dir, save_date, save_time, save_frequency, save_height, save_meanpower):
    """
    Converts MD2/MD4 radar data files to CSV format, extracting key features.

    Args:
        input_file (str): Path to the input MD2/MD4 file.
        output_dir (str): Path to the directory where the converted CSV file will be saved.
        save_date (bool): Flag to save date information.
        save_time (bool): Flag to save time information.
        save_frequency (bool): Flag to save frequency information.
        save_height (bool): Flag to save height information.
        save_meanpower (bool): Flag to save mean power information.

    Returns:
        None
    """
    max_ntimes = 256
    max_ndopbins = 300000
    dheight = 3.0  # not defined in data file, from documentation p.

    with open(input_file, "rb") as f:
        f.seek(-1, 2)  # go to the file end.
        eof = f.tell()  # get the end of file location
        f.seek(0, 0)  # go back to file beginning

        # 1) read header information as described in the documentation p. 26-27
        site = f.read(3).decode("utf-8")
        ascii_datetime = f.read(22).decode("utf-8")
        filetype = f.read(1).decode("utf-8")
        nfreqs = struct.unpack("<H", f.read(2))[0]
        ndops = struct.unpack("<B", f.read(1))[0]
        minheight = struct.unpack("<H", f.read(2))[0]
        maxheight = struct.unpack("<H", f.read(2))[0]
        pps = struct.unpack("<B", f.read(1))[0]
        npulses_avgd = struct.unpack("<B", f.read(1))[0]
        base_thr100 = struct.unpack("<H", f.read(2))[0]
        noise_thr100 = struct.unpack("<H", f.read(2))[0]
        min_dop_forsave = struct.unpack("<B", f.read(1))[0]
        dtime = struct.unpack("<H", f.read(2))[0]
        gain_control = f.read(1).decode("utf-8")
        sig_process = f.read(1).decode("utf-8")
        noofreceivers = struct.unpack("<B", f.read(1))[0]
        spares = f.read(11).decode("utf-8")

        month = ascii_datetime[1:4]
        day = int(ascii_datetime[5:7])
        hour = int(ascii_datetime[8:10])
        minute = int(ascii_datetime[11:13])
        sec = int(ascii_datetime[14:16])
        year = int(ascii_datetime[17:21])

        month_number = datetime.datetime.strptime(month, '%b').month
        mydate = datetime.date(year, month_number, day)
        jd = mydate.toordinal() + 1721424.5
        jd0jd = datetime.date(1986, 1, 1)
        jd0 = jd0jd.toordinal() + 1721424.5

        time_header = (jd - jd0) * 86400 + hour * 3600 + minute * 60 + sec
        time_hour = 3600 * (time_header / 3600)

        # 2) read all frequencies used
        freqs = struct.unpack("<" + "f" * nfreqs, f.read(4 * nfreqs))

        if filetype == 'I':
            max_nfrebins = nfreqs
        else:
            max_nfrebins = min(max_ntimes * nfreqs, max_ndopbins)

        # 3) read rawdata
        nheights = int(maxheight / dheight + 1)

        times = []
        frebins = []
        frebins_x = []
        frebins_gain_flag = []
        frebins_noise_flag = []
        frebins_noise_power10 = []
        time_min = 0
        time_sec = 0
        timex = -1
        freqx = nfreqs - 1
        dopbinx = -1
        frebinx = -1
        iq_bytes = np.zeros((noofreceivers, 2))
        dopbin_x_timex = []
        dopbin_x_freqx = []
        dopbin_x_hflag = []
        dopbin_x_dop_flag = []
        dopbin_iq = []
        hflag = 0

        time_min = struct.unpack("<B", f.read(1))[0]
        while time_min != 255:
            time_sec = struct.unpack("<B", f.read(1))[0]
            flag = struct.unpack("<B", f.read(1))[0]  # gainflag
            timex += 1
            times.append(time_hour + 60 * time_min + time_sec)
            for freqx in range(nfreqs):
                noise_flag = struct.unpack("<B", f.read(1))[0]  # noiseflag
                noise_power10 = struct.unpack("<H", f.read(2))[0]
                frebinx += 1
                frebins_gain_flag.append(flag)
                frebins_noise_flag.append(noise_flag)
                frebins_noise_power10.append(noise_power10)
                flag = struct.unpack("<B", f.read(1))[0]
                while flag < 224:
                    ndops_oneh = struct.unpack("<B", f.read(1))[0]
                    hflag = flag
                    if ndops_oneh >= 128:
                        ndops_oneh = ndops_oneh - 128
                        hflag = hflag + 200
                    for dopx in range(ndops_oneh):
                        dop_flag = struct.unpack("<B", f.read(1))[0]
                        for rec in range(noofreceivers):
                            iq_bytes[rec, 0] = struct.unpack("<B", f.read(1))[0]
                            iq_bytes[rec, 1] = struct.unpack("<B", f.read(1))[0]
                        dopbinx += 1
                        dopbin_iq.append(iq_bytes.copy())
                        dopbin_x_timex.append(timex)
                        dopbin_x_freqx.append(freqx)
                        dopbin_x_hflag.append(hflag)
                        dopbin_x_dop_flag.append(dop_flag)
                    flag = struct.unpack("<B", f.read(1))[0]  # next hflag/gainflag/FF
            time_min = flag
            if (f.tell() - 1) != eof:
                time_min = struct.unpack("<B", f.read(1))[0]  # next record

    meanpower = []
    frequency = []
    # calculate average power
    for idx in range(len(dopbin_iq)):
        # transform coordinates
        for rec in range(noofreceivers):
            for comp in range(2):
                if dopbin_iq[idx][rec][comp] > 127:
                    dopbin_iq[idx][rec][comp] = dopbin_iq[idx][rec][comp] - 256

        absvalue1 = math.sqrt((dopbin_iq[idx][0][0]) ** 2 + (dopbin_iq[idx][0][1]) ** 2)
        absvalue2 = math.sqrt((dopbin_iq[idx][1][0]) ** 2 + (dopbin_iq[idx][1][1]) ** 2)
        absvalue3 = math.sqrt((dopbin_iq[idx][2][0]) ** 2 + (dopbin_iq[idx][2][1]) ** 2)
        absvalue4 = math.sqrt((dopbin_iq[idx][3][0]) ** 2 + (dopbin_iq[idx][3][1]) ** 2)
        mvalue = median([absvalue1, absvalue2, absvalue3, absvalue4])
        # mvalue = (absvalue1 + absvalue2 + absvalue3 + absvalue4) / 4.
        if mvalue == 0:
            power = 0
        else:
            power = 20 * math.log10(mvalue)
        meanpower.append(power)
        frequency.append(freqs[dopbin_x_freqx[idx]] / 1000000.0)

    height = list(np.array(dopbin_x_hflag) * 3)

    date = [str(day) + "/" + str(month_number).zfill(2) + "/" + str(year)] * len(frequency)
    time = [ascii_datetime[8:16]] * len(frequency)

    # Initialize lists for selected features
    selected_features = []
    if save_date:
        selected_features.append(date)
        header = ["Date"]
    if save_time:
        selected_features.append(time)
        header = header + ["Time"]
    if save_frequency:
        selected_features.append(frequency)
        header = header + ["Frequency"]
    if save_height:
        selected_features.append(height)
        header = header + ["Height"]
    if save_meanpower:
        selected_features.append(meanpower)
        header = header + ["MeanPower"]

    # Transpose the selected features to have the same length
    selected_features = list(map(list, zip(*selected_features)))

    # Create a CSV file and write the selected features
    csv_filename = os.path.splitext(os.path.basename(input_file))[0] + ".csv"
    csv_path = os.path.join(output_dir, csv_filename)
    with open(csv_path, mode="w", newline="") as csv_file:
        csv_writer = csv.writer(csv_file)
        if selected_features:
            csv_writer.writerow(header)
            for row in selected_features:
                csv_writer.writerow(row)


# Function to log information about corrupt files
def log_corrupt_file(file_path, error_message):
    """
    Logs information about a corrupt input file.

    Args:
        file_path (str): The path to the corrupt input file.
        error_message (str): The error message describing the issue.

    Returns:
        None
    """
    with open("corrupt_files.txt", "a") as log_file:
        log_file.write(f"{file_path}: {error_message}\n")


# Function to find MD2 and MD4 files within a folder and its subdirectories
def find_md2_md4_files(folder):
    """
    Finds MD2 and MD4 files within a folder and its subdirectories.

    Args:
        folder (str): The path to the folder to search for MD2/MD4 files.

    Returns:
        list: A list of paths to MD2/MD4 files found.
    """
    md2_md4_files = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.lower().endswith((".md2", ".md4")):
                md2_md4_files.append(os.path.join(root, file))
    return md2_md4_files


# Main function that handles command-line arguments and initiates the conversion
def main():
    parser = argparse.ArgumentParser(description="MDX2CSV Converter - Extracts MD2/MD4 key features to CSV",
                                     usage="mdx2csv_cli.py [command args]",
                                     prog="MDX2CSV_version 1")
    parser.add_argument("input_path", type=str, help="Input folder or file to convert")
    parser.add_argument("output_directory", type=str, help="Output directory for converted files")
    parser.add_argument("--custom", type=str, default="dtfhm",
                        help="Specify custom features to save using a shorthand notation (e.g., 'dtfhm' for Date, Time,"
                             " Frequency, Height, MeanPower)")
    args = parser.parse_args()

    input_path = args.input_path  # Use either positional or optional argument
    output_directory = args.output_directory  # Use either positional or optional argument

    if not output_directory:
        print("Specify an output directory.")
        return

    if os.path.isdir(input_path):
        input_files = find_md2_md4_files(input_path)
    else:
        input_files = [input_path]

    if not input_files:
        print("No MD2 or MD4 files found in the specified folder.")
        return

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Define feature selection based on --custom argument
    custom_features = set(args.custom)

    for input_file in tqdm(input_files):
        try:
            mdx2csv(input_file, output_directory, 'd' in custom_features, 't' in custom_features,
                    'f' in custom_features, 'h' in custom_features, 'm' in custom_features)
        except Exception as e:
            log_corrupt_file(input_file, str(e))

    print("Conversion completed!")


if __name__ == "__main__":
    main()
