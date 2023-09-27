import argparse
import os
import struct
import datetime
import numpy as np
from tqdm import tqdm


# Function to perform the actual MD2/MD4 to TXT conversion
def mdx2txt(input_file, output_dir):
    """
    Converts an MD2/MD4 input file to TXT format.

    Args:
        input_file (str): The path to the input MD2/MD4 file.
        output_dir (str): The path to the directory where the converted file will be saved.

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

    # Open the output file in write mode
    fname = os.path.join(output_dir, input_file.split("\\")[-1][:-4] + ".txt")
    with open(fname, "w") as output_f:
        # Writing the header information to the output file
        output_f.write(f"Site: {site}\n")
        output_f.write(f"Date/Time: {ascii_datetime}\n")
        output_f.write(f"Filetype: {filetype}\n")
        output_f.write(f"Number of Frequencies: {nfreqs}\n")
        output_f.write(f"Number of Doppler Bins: {ndops}\n")
        output_f.write(f"Minimum Height: {minheight}\n")
        output_f.write(f"Maximum Height: {maxheight}\n")
        output_f.write(f"Height Interval: {dheight}\n")  # not in data file
        output_f.write(f"Pulses per Second: {pps}\n")
        output_f.write(f"Number of Pulses Averaged: {npulses_avgd}\n")
        output_f.write(f"Base Threshold (100): {base_thr100}\n")
        output_f.write(f"Noise Threshold (100): {noise_thr100}\n")
        output_f.write(f"Minimum Doppler for Save: {min_dop_forsave}\n")
        output_f.write(f"Time Between Measurements (seconds): {dtime}\n")
        output_f.write(f"Gain Control: {gain_control}\n")
        output_f.write(f"Signal Processing: {sig_process}\n")
        output_f.write(f"Number of Receivers: {noofreceivers}\n")
        output_f.write(f" mdx2txt Version: 1.0 #COMMENT\n")
        # Writing all frequencies used
        output_f.write("Frequencies Used:\n")
        for freq in freqs:
            output_f.write(f"{freq}\n")

        # Writing the raw data to the output file
        output_f.write("Raw Data:\n")
        for i in range(len(times)):
            output_f.write(f"Time: {times[i]}, Min: {time_min}, Sec: {time_sec}\n")
            for freqx in range(nfreqs):
                # Writing data for each frequency
                output_f.write(f"Frequency: {freqx}, Gain Flag: {frebins_gain_flag[freqx]}, "
                               f"Noise Flag: {frebins_noise_flag[freqx]}, "
                               f"Noise Power (10x): {frebins_noise_power10[freqx]}\n")
                for dopx in range(ndops):
                    # Writing data for each Doppler bin
                    dopbin_index = dopx + freqx * ndops
                    iq_values = dopbin_iq[dopbin_index][1]
                    formatted_iq = f"[{iq_values[0]:<5}, {iq_values[1]:<5}]"
                    output_f.write(f"Doppler Bin: {dopx}, Receiver IQ: {formatted_iq}\n")


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
    parser = argparse.ArgumentParser(description="MDX2TXT Converter - Convert MD2/MD4 files to TXT format",
                                     usage="mdx2txt_cli.py [command args]",
                                     prog="MDX2ASCII_version 1")
    parser.add_argument("input_path", type=str, help="Input folder or file to convert")
    parser.add_argument("output_directory", type=str, help="Output directory for converted files")
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

    for input_file in tqdm(input_files):
        try:
            mdx2txt(input_file, output_directory)
        except Exception as e:
            log_corrupt_file(input_file, str(e))

    print("Conversion completed!")


if __name__ == "__main__":
    main()
