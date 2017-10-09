__author__ = 'iob'

if __name__ == "__main__":
    __author__ = 'Srikanth'
    parser = argparse.ArgumentParser(description='Program to parse...')
    parser.add_argument('-d', '--database', help='Input db', required=True)
    parser.add_argument('-i', '--inputfile', help='Input file', required=True)
    #parser.add_argument('-o','--output',help='Output file name', required=True)
    args = parser.parse_args()