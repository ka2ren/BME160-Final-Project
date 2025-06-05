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


def get_positive_strand(GRCH38_file, max_records=5):
    ''' Extract sequences and ranges for positive strand from GRCh38 FASTA file. '''
    pos_strand_seqs = []
    pos_strand_ranges = []
    reader = FastAreader(GRCH38_file)
    count = 0
    for header, sequence in reader.readFasta():
        if "strand=+" in header:
            # Extract range from header
            parts = header.split()
            for part in parts:
                if part.startswith("range="):
                    range_str = part.split("=")[1]
                    chrom, coords = range_str.split(":")
                    start, end = map(int, coords.split("-"))
                    pos_strand_ranges.append((start, end))
            pos_strand_seqs.append(sequence)
            count += 1
            if count >= max_records:
                break
    return pos_strand_seqs, pos_strand_ranges


def extract_t2t_regions(t2t_fasta, ranges):
    # Assumes T2T fasta has a single full chr1 sequence
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
    """Count the number of telomeric repeats in a sequence."""
    return seq.upper().count(repeat)

def compare_telomeric_repeats(grch38_seqs, t2t_regions, repeat="TTAGGG"):
    print(f"{'Index':<5} {'GRCh38_count':<15} {'T2T_count':<10}")
    print("-" * 35)
    for i, (gseq, tseq) in enumerate(zip(grch38_seqs, t2t_regions)):
        g_count = count_telomeric_repeats(gseq, repeat)
        t_count = count_telomeric_repeats(tseq, repeat)
        print(f"{i:<5} {g_count:<15} {t_count:<10}")

if __name__ == "__main__":
    grch38_seqs, grch38_ranges = get_positive_strand("GRCH38_test.fasta", max_records=5)
    t2t_regions = extract_t2t_regions("T2T_test.fasta", grch38_ranges)
    compare_telomeric_repeats(grch38_seqs, t2t_regions)



 
    