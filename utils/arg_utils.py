import argparse


def add_common_params(parser: argparse.ArgumentParser, default_outfile: str):
    parser.add_argument(
        "end",
        type=str,
        nargs="?",
        default="120",
        help="If a number: represents the number of days after the start date to calculate for, otherwise attempt to interpert as a date with format %%Y-%%m-%%d; default of 120 days",
    )

    parser.add_argument(
        "-tz",
        "--timezone",
        type=str,
        default="UTC",
        help="The timezone to calculate with respect to, default is UTC",
    )

    parser.add_argument(
        "-r",
        "--reload",
        action="store_true",
        help="Re-download the models file from the internet",
    )

    parser.add_argument(
        "-ca",
        "--cache",
        action="store_true",
        help="Automatically Cache model files with the format of %%Y%%m%%d_sats.txt",
    )

    parser.add_argument(
        "-l",
        "--all",
        action="store_true",
        help="Use all satellites from database; recommended to use if --tles is also used",
    )

    parser.add_argument(
        "-sd",
        "--start-date",
        type=str,
        default="today",
        help="specify a starting date to begin calculations; format is %%Y-%%m-%%d",
    )

    parser.add_argument(
        "-il",
        "--ignore-limit",
        action="store_true",
        help="ignore the 14 day recommended TLE limit on satellite data",
    )

    parser.add_argument(
        "--tles",
        type=str,
        nargs="+",
        help="Specify the model files to use instead of downloading them; supports TLE, OMM CSV, and OMM XML formats",
    )

    parser.add_argument(
        "--chain",
        action="store_true",
        help="When the inputted OMM or TLE file has many models in sequence for various satellites - chain them together",
    )

    parser.add_argument(
        "-o",
        "--output-file",
        type=str,
        help="name of the file to output",
        default=default_outfile,
    )
