import sys
import argparse
from sage.all import *
from Kempe_sep_class import Configuration

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()
    config = Configuration(6, args.filename)
    config.is_reducible()

if __name__ == '__main__':
    main()