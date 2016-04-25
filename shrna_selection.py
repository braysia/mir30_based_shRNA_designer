from __future__ import division
from Bio.Seq import Seq
import bioservices as bs
import re
from splinter import Browser
import time
from Bio import Entrez


email = '910taka@gmail.com'
sub_email = ''

class CdsGetter(object):
    def __init__(self, nm_id):
        self.nm_id = nm_id
        # self.eut = bs.EUtils()

    def run(self):
        self.retrieve_gb_entrez()
        gb_sequences = self.retrieve_sequences()
        start_pos, end_pos = self.retrieve_coding_pos()

        print gb_sequences, start_pos, end_pos
        cds_seq = gb_sequences[start_pos-1:end_pos]
        return cds_seq

    def retrieve_sequences(self):
        # fasta_seq = self.eut.EFetch('nucleotide', self.nm_id, rettype='fasta')
        # sequences = fasta_seq.split('\n')[1:]
        sequences = self.extract_origin(self.gb_seq)
        return sequences

    def retrieve_coding_pos(self):
        # gb_seq = self.eut.EFetch('nucleotide', self.nm_id, rettype='gb')
        gb_seq = self.gb_seq
        gb_seq = gb_seq.split('\n')
        cds_line = [i for i in gb_seq if '     CDS' in i]
        assert len(cds_line) == 1
        regexp = re.search('[A-Z, a-z]*(?P<fir>[0-9]*)..(?P<sec>[0-9]*)', cds_line[0])
        start_pos, end_pos = int(regexp.group('fir')), int(regexp.group('sec'))
        return start_pos, end_pos

    def retrieve_gb_entrez(self):
        Entrez.email = email
        handle = Entrez.efetch('nucleotide', id=self.nm_id, rettype='gb', retmode='text')
        gb_seq = handle.readlines()
        gb_seq = ''.join(gb_seq)
        self.gb_seq = gb_seq

    @staticmethod
    def extract_origin(gb_seq):
        '''return ORIGIN sequences'''
        gb = gb_seq.split('ORIGIN')[-1]
        seq = []
        for i in gb.split('\n'):
            seq.append(i[10:].replace(" ", ""))
        seq = ''.join(seq)
        return seq


def pick_rna(seq):
    '''Put anti-sense (anti-sense) sequence as an input.
    Based on Nature Protocol 2012 Lucas Dow...Scott Lowe.
    '''
    assert len(seq) == 21

    print '1. 1st position is {0}'.format(seq[0])
    assert seq[0]=='A' or seq[0]=='U'
    au_ratio = len([i for i in seq if i=='A' or i=='U'])/len(seq)
    print '2. A/U content is {0}'.format(au_ratio)
    assert 0.4 <= au_ratio <= 0.8

    au_ratio1to14 = len([i for i in seq[:14] if i=='A' or i=='U'])/len(seq[:14])
    print '3. A/U content from 1 through 14 is {0}'.format(au_ratio1to14)
    assert 0.5 <= au_ratio1to14
    au_ratio15to21 = len([i for i in seq[14:] if i=='A' or i=='U'])/len(seq[14:])
    print '4. AU1-14/AU15-21 is {0}'.format(au_ratio1to14/au_ratio15to21)
    assert 1 <= au_ratio1to14/au_ratio15to21
    print '5. 20th position is {0}'.format(seq[19])
    assert seq[20] != 'A'
    isau = seq[12]=='A' or seq[12] == 'U'
    isu = seq[13] == 'U'
    print '6. 13th and 14th positions are {0} and {1}'.format(seq[12], seq[13])
    assert isau or isu
    for i in ['AAAAAA', 'TTTTT', 'CCCC', 'GGGG']:
        assert seq.find(i) == -1
    print '7. not containing AAAAAA TTTTT CCCC GGGG'
    print '\nguide strand: {0}'.format(seq)
    print 'revcon strand: {0}'.format(Seq(seq).reverse_complement())

class AutomatedSirnaDesigner(object):
    """docstring for ClassName"""
    def __init__(self, nm_id, cds_seq):
        self.nm_id = nm_id
        self.cds_seq = cds_seq

    def run(self):
        html = self.run_analysis()
        predicted = self.extract_html(html)
        return predicted

    def run_analysis(self):
        url = 'http://biodev.extra.cea.fr/DSIR/DSIR.html'
 
        with Browser() as browser:
            browser.visit(url)
            browser.find_by_name('seqname').fill('test')
            browser.find_by_name('sequence_all').fill('>test\n'+self.cds_seq)
            browser.find_by_value('Run analysis').click()
            result = browser.html
        return result

    def extract_html(self, html):
        splitted = html.split('\n')
        table = [i for i in splitted if 'input type="checkbox"' in i][0]
        table = table.split('align="center')[2:]
        pattern = 'value."(?P<value>[0-9]*)"[^0-9]*[0-9]*[^0-9]*(?P<pos>[0-9]*)[^A-Z]*(?P<sense>[A-Z]*)[^A-Z]*(?P<antisense>[A-Z]*)'
        store = []
        for t in table:
            regexp = re.search(pattern, t)
            store.append((regexp.group('value'), regexp.group('pos'), regexp.group('sense'), regexp.group('antisense')))
        return store            


class AddFifthOftheStrand(object):
    '''predicted_set is a tuple of number, position, sense and antisense sequences without 5'th of the sense
    and 3' of the guide strand. This class will translate 21-bp to 22-bp. 
    '''
    def __init__(self, cds_seq, predicted_set):
        self.cds_seq = cds_seq
        self.predicted_set = [list(i) for i in predicted_set] # turn tuple to list

    def run(self):
        for pset in self.predicted_set:
            self.turn_21to22(self.cds_seq, pset)
        return self.predicted_set

    def turn_21to22(self, cds_seq, pset):
        pos = int(pset[1])
        base = Seq(cds_seq).transcribe()[pos-2]
        if base in ['A', 'U', 'a', 'u']:
            pset[2] = 'C' + pset[2]
            pset[3] = pset[3] + str(Seq(base).reverse_complement())
        elif base in ['C', 'G', 'c', 'g']:
            pset[2] = 'A' + pset[2]
            pset[3] = pset[3] + str(Seq(base).reverse_complement())


def turn_u_to_t(seq):
    seq.replace('U', 'T')
    seq.replace('u', 'T')
    return seq

if __name__ == '__main__':
    # seq = 'UUGCGUAGCACAAAUUUCGGU'
    # pick_rna(seq)

    # nm_id = 'NM_007912.4'
    # cds_seq = CdsGetter(nm_id).run()
    # Dsir(nm_id, cds_seq).run()
    pass