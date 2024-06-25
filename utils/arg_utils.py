import argparse


def add_common_params(parser: argparse.ArgumentParser, default_outfile: str):
    parser.add_argument(
        "days",
        type=int,
        nargs="?",
        default=120,
        help="The number of days to calculate for, default of 120 days",
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
        "-ca", "--cache", action="store_true", help="Automatically Cache model files with the format of %%Y%%m%%d_sats.txt"
    )

    parser.add_argument(
        "-l", "--all", action="store_true", help="Use all satellites from database"
    )

    parser.add_argument(
        "-sd",
        "--start-date",
        type=str,
        default="today",
        help="specify a starting date to begin calculations format is <year>-<month>-<day>",
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
        "-o",
        "--output-file",
        type=str,
        help="name of the file to output",
        default=default_outfile,
    )
