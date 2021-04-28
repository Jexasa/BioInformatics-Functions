from BioStructs import DNA_Codons, Nucleotides
import random
from collections import Counter

class Bio_sequence : 

    def __init__(self,seq="ATCG",seqtype="DNA",label="NoLabel"):
        """Initiliazation + validation  """
        self.seq= seq.upper()
        self.seqtype = seqtype
        self.label = label
        self.is_valid = self.__Validate()
        assert self.is_valid, f" incorrect data provided {self.seqtype} sequence"

    def __Validate(self):
        """ Validates the sequence """

        #We check if the nucleotides in the sequence are in our nucleotide list
        return set(Nucleotides).issuperset(self.seq)

    def Get_sequence_type(self):
        """ returns the type of sequence"""
        return self.seqtype

    def Information(self):
        """ Full sequence information """
        return f"-Label- {self.label}\n-Sequence- {self.seq}\n-Type- {self.seqtype}\n-Length-{len(self.seq)}\n"

    def Generate_random_dna_sequence(self,length=42,seqtype="DNA"):
        """Generate a random dna sequence and reinitializes it """
        seq = ''.join([random.choice(Nucleotides) for i in range (length)])
        self.__init__(seq,seqtype,"Random Sequence")

    def Nucleotides_frequency(self):
        """ Returns the number of nucleotides in a dna sequence """
        return dict(Counter(self.seq))

    def Transcript(self):
        """ Replaces Thymine nucleotides with Uracile nucleotides
            aka Transcription"""
        return self.seq.replace("T","U")

    def ReverseComplement(self):
        """ Swaps A with T and G with C.
        Then reverses and returns the string"""

        #mapping contains the changed sequence
        mapping = str.maketrans('ATCG','TAGC')

        #reverse the string of the mutated sequence
        return self.seq.translate(mapping)[::-1]

    
    def GCContent(self):
        """ Returns the appearance frequency of G and C of the sequence""" 
        g = 0
        c = 0
        for i in self.seq:
            if  i == 'G':
                g+=1
            elif i == 'C':
                c+=1
        freq= g + c / len(self.seq) * 100
        rounded_freq = (f'{freq:.3f}')
        return rounded_freq

    def CContentSubsec(self,k=20):
        """Returns the appearance frequency of G and C of the subsequences
            k=the length of the subsequence """
        subsecs=[]
        
        for i in range (0, len(self.seq) - k + 1, k):
            
            #Cuts the list from i to i+k
            subseq = self.seq[i:i +k]

            #appends the frequency of the subseq
            subsecs.append(
                round(subseq.count('C') + subseq.count('G') / len(subseq) * 100))
        return subsecs

    def TranslateSequence(self,init_pos=0):
        """returns the aminoacids of the sequence"""

        return [DNA_Codons[self.seq[pos:pos+3]] for pos in range (init_pos,len(self.seq) - 2,3)]

    def ReadFrames (self):

        """Creates a list of all (6) open reading frames of a given sequence """
        frames = []
        for i in range (3):
            frames.append(self.TranslateSequence(i))

        tmp = Bio_sequence(self.ReverseComplement(), self.seqtype)
        for i in range (3):
            frames.append(tmp.TranslateSequence(i))
        return frames

    def Proteins_From_ReadFrame(self,sequence):
        """ Read proteins from a readframe sequence and returns a list from start "M" to "_" """
        tmp = []
        proteins=[]

        #for every aminoacid-codon in the given sequence
        for aminoacid in sequence :
            if aminoacid == "_" :
                # tmp !epmty
                if tmp :
                    #for every start codon
                    for i in tmp :
                        proteins.append(i)
            else :
                if aminoacid == "M":
                    #We can have to M's so we seperate them
                    tmp.append("")
                for i in range (len(tmp)):
                    tmp[i] += aminoacid
        return proteins


    def colored(self):

        """ Colors the nucleotides of a dna sequence """
        colors = {
            'A': '\033[92m',
            'C': '\033[94m',
            'G': '\033[93m',
            'T': '\033[91m',
            'U': '\033[91m',
            'reset': '\033[0;0m' }
        tmpstr =''

        for n in self.seq:
            if n in colors :
                tmpstr += colors[n] + n 
            else :
                tmpstr += colors['reset'] + n
        
        return tmpstr + '\033[0;0m'