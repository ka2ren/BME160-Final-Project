import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence


def get_positive_strand(GRCH38_file):
    ''' Extract sequences and ranges for positive strand from GRCh38 FASTA file. '''
    pos_strand_seqs = []
    pos_strand_ranges = []
    reader = FastAreader(GRCH38_file)
    
    for header, sequence in reader.readFasta():
        if "strand=+" in header:
            parts = header.split()
            for part in parts:
                if part.startswith("range="):
                    range_str = part.split("=")[1]
                    chrom, coords = range_str.split(":")
                    start, end = map(int, coords.split("-"))
                    pos_strand_ranges.append((start, end))
            pos_strand_seqs.append(sequence)

    return pos_strand_seqs, pos_strand_ranges


def extract_t2t_regions(t2t_fasta, ranges):
    # T2T fasta has a single full chr1 sequence
    reader = FastAreader(t2t_fasta)
    for header, sequence in reader.readFasta():
        t2t_seq = sequence
        break  # Only one sequence expected
    matched_regions = []
    for start, end in ranges:
        # FASTA is 1-based, Python is 0-based
        matched_regions.append(t2t_seq[start-1:end])
    return matched_regions


def count_telomeric_repeats(seq, repeat="TTAGGG"):
    """Count overlapping occurrences of telomeric repeats in a sequence."""
    seq = seq.upper()
    repeat = repeat.upper()
    count = 0
    i = 0
    while i <= len(seq) - len(repeat):
        if seq[i:i+len(repeat)] == repeat:
            count += 1
            i += 1  # Move by 1 to allow overlaps
        else:
            i += 1
    return count

def count_terminal_telomeric_repeats(fasta_file, region_size=10000, repeat="TTAGGG"):
    """Returns (header, five_prime_count, three_prime_count) for the first sequence in the file."""
    reader = FastAreader(fasta_file)
    for header, sequence in reader.readFasta():
        seq = sequence.upper()
        five_prime = seq[:region_size]
        three_prime = seq[-region_size:]
        five_count = count_telomeric_repeats(five_prime, repeat)
        three_count = count_telomeric_repeats(three_prime, repeat)
        return header, five_count, three_count  # Only process the first sequence

if __name__ == "__main__":
    t2t_header, t2t_5, t2t_3 = count_terminal_telomeric_repeats("T2T_test.fasta", region_size=10000, repeat="TTAGGG")
    grch38_header, grch38_5, grch38_3 = count_terminal_telomeric_repeats("GRCH38_test.fasta", region_size=10000, repeat="TTAGGG")
    print(f"T2T: counts 5': {t2t_5}; 3': {t2t_3}   \nGRCH38: counts 5': {grch38_5}; 3': {grch38_3}")

    # Calculate calibration factors (ratios)
    cal_factor_5 = t2t_5 / grch38_5 if grch38_5 != 0 else float('inf')
    cal_factor_3 = t2t_3 / grch38_3 if grch38_3 != 0 else float('inf')
    print(f"Calibration factor (T2T/GRCH38) 5': {cal_factor_5:.2f}, 3': {cal_factor_3:.2f}")

'''
Output:
range: 248368735-248511794 T2T: 420 GRCH38: 24
'''

 
    