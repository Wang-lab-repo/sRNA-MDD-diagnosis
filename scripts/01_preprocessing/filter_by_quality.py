import gzip
import sys

def read_fastq(fh):
    while True:
        try:
            header = next(fh).strip()
            seq = next(fh).strip()
            plus = next(fh).strip()
            qual = next(fh).strip()
            yield header, seq, plus, qual
        except StopIteration:
            break

def count_low_qual(qual_str, threshold):
    return sum([ord(c) - 33 < threshold for c in qual_str])

def filter_reads(input_file, output_file):
    opener = gzip.open if input_file.endswith('.gz') else open
    with opener(input_file, 'rt') as in_fh, open(output_file, 'w') as out_fh:
        for header, seq, plus, qual in read_fastq(in_fh):
            q10_count = count_low_qual(qual, 10)
            q13_count = count_low_qual(qual, 13)
            if q10_count > 4 or q13_count > 6:
                continue
            out_fh.write(f"{header}\n{seq}\n{plus}\n{qual}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python filter_by_quality.py input.fq[.gz] output.fq")
        sys.exit(1)
    filter_reads(sys.argv[1], sys.argv[2])
