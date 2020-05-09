from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

class PCRSimulator:
    def __init__(self, for_primer, rev_primer, template):
        self.for_primer = for_primer
        self.rev_primer = rev_primer
        self.template = template

    DNA_bases = ['a', 't', 'c','g']
    DNA_degen = ['a', 't', 'c','g', 'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n']
    rev_comp_dic = {'a': 't',
                    't': 'a',
                    'g': 'c',
                    'c': 'g',
                    'r': 'y',
                    'y': 'r',
                    's': 's',
                    'w': 'w',
                    'k': 'm',
                    'm': 'k',
                    'n': 'n'}

    # Check to see if string has only DNA bases.
    def only_dna(seq):
        for base in seq:
            if base not in PCRSimulator.DNA_bases:
                return False
        return True

    # Check to see if string has only DNA bases or degeneracy codes.
    def only_dna_or_degen(seq):
        for base in seq:
            if base not in PCRSimulator.DNA_degen:
                return False
        return True

    # Reverse Complement Calculator.
    def rev_comp_calc(seq):
        output = ''
        for base in seq:
            output += PCRSimulator.rev_comp_dic[base]
        output = output[::-1]
        return output

    # Finding the annealing region for the forward primer. Raises error if it does not anneal.
    def find_anneal_reg_for(oligo, seq):
        x = 6 # Minimum length required for binding.
        anneal_seq, temp = oligo[len(oligo)-x:], oligo[len(oligo)-x:]
        # Check to see if it is the forward primer.
        if anneal_seq in seq:
            while temp in seq and x <= len(oligo):
                anneal_seq = temp
                anneal_start = seq.find(anneal_seq)
                x += 1
                if x > len(oligo):
                    return [anneal_start, anneal_seq]
                temp = oligo[len(oligo)-x:]
            # Check to see if there are multiple instances of primer.
        else:
            raise Exception('Primer does not have homologous sequence to template strain.')
        PCRSimulator.check_annealing_strength(anneal_seq)
        return [anneal_start, anneal_seq]

    def find_anneal_reg_rev(oligo, seq):
        x = 6
        oligo = PCRSimulator.rev_comp_calc(oligo)
        anneal_seq, temp = oligo[:6], oligo[:6]
        if anneal_seq in seq:
            while temp in seq and x <= len(oligo):
                anneal_seq = temp
                anneal_start = seq.find(anneal_seq)
                x += 1
                if x > len(oligo):
                    return [anneal_start, anneal_seq]
                temp = oligo[:x]
            # Check to see if there are multiple instances of primer.
        else:
            raise Exception('Primer does not have homologous sequence to template strain.')
        PCRSimulator.check_annealing_strength(anneal_seq)
        return [anneal_start, PCRSimulator.rev_comp_calc(anneal_seq)]

    # Resets the start of the plasmid to a certain base position.
    def circularize(seq, start):
        seq = seq[start:] + seq[:start]
        return seq

    def trim_edges(self):
        for_anneal = PCRSimulator.find_anneal_reg_for(self.for_primer, self.template)
        self.template = PCRSimulator.circularize(self.template, for_anneal[0])
        self.template = self.for_primer + self.template[len(for_anneal[1]):]
        rev_anneal = PCRSimulator.find_anneal_reg_rev(self.rev_primer, self.template)
        template = self.template[:rev_anneal[0]] + PCRSimulator.rev_comp_calc(self.rev_primer)
        return template

    def pcr_simulator(self):
        self.fix_oligos()
        self.check_primers()
        product = self.trim_edges()
        return product

    def fix_oligos(self):
        self.for_primer = self.for_primer[:len(self.for_primer) - 1].lower()
        self.rev_primer = self.rev_primer[:len(self.rev_primer) - 1].lower()
        self.template = self.template[:len(self.template) - 1].lower()
        self.for_primer = self.for_primer.replace('\n', '')
        self.rev_primer = self.rev_primer.replace('\n', '')
        self.template = self.template.replace('\n', '')
        self.for_primer = self.for_primer.replace(' ', '')
        self.rev_primer = self.rev_primer.replace(' ', '')
        self.template = self.template.replace(' ', '')

    # Test to see if the primer's melting temperature is high enough to undergo PCR.
    # Raises error if primer is not good enough.
    def bad_primer(oligo):
        if mt.Tm_NN(Seq(oligo)) < 50:
            raise Exception('''Primer has too low of melting temperature to work in
            PCR; try increasing length or GC content.''')

    def check_primers(self):
        if not PCRSimulator.only_dna_or_degen(self.for_primer) or not PCRSimulator.only_dna_or_degen(self.rev_primer) or not PCRSimulator.only_dna_or_degen(self.template):
            raise Exception('Invalid characters detected in input.')
        else:
            PCRSimulator.bad_primer(self.for_primer)
            PCRSimulator.bad_primer(self.rev_primer)

    def check_annealing_strength(seq):
        if not PCRSimulator.only_dna(seq):
            raise Exception('Degeneracy codes detected in annealing region of primer.')

