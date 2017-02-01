from __future__ import print_function

import six
import copy
import os
import unittest
import sys
import pandas as pd
from pandas.util.testing import assert_frame_equal

from tracerlib import base_dir
from tracerlib.io import parse_IgBLAST, parse_invariant_cells
from tracerlib.tasks import Assembler


class TestAssemble(unittest.TestCase):

    expected_folder = os.path.join(base_dir, 'test_data', 'expected_summary')
    results_folder = os.path.join(base_dir, 'test_data', 'results',
                                  'filtered_TCR_summary')

    def setUp(self):
        config_file = os.path.expanduser('~/.tracerrc')
        assert os.path.isfile(config_file), "Config file ~/.tracerrc not found"
        fastq1 = os.path.join(base_dir, 'test_data', 'cell1_1.fastq')
        fastq2 = os.path.join(base_dir, 'test_data', 'cell1_2.fastq')
        self.assemble_args = dict(
            ncores='1', config_file=config_file,
            resume_with_existing_files=False, species='Mmus',
            seq_method='imgt', fastq1=fastq1, fastq2=fastq2,
            cell_name='cell1', output_dir=self.results_folder, single_end=False,
            fragment_length=False, fragment_sd=False, quant_method='salmon')
        self.assembler = Assembler(**self.assemble_args)

    def test_args(self):

        assemble_args = copy.deepcopy(self.assemble_args)
        assemble_args["fastq1"] = '/not/a/file'
        assemble_args["fastq2"] = '/not/a/file'
        # Launcher checks for fastq file existance
        with self.assertRaisesRegexp(OSError, 'FASTQ file not found'):
            Assembler(**assemble_args).run()

        # Launcher checks for second fastq if paired end
        assemble_args['fastq2'] = None
        assemble_args['fastq1'] = os.path.join(base_dir, 'test_data',
                                               'cell1_1.fastq')
        with self.assertRaisesRegexp(AssertionError, 'Only one fastq'):
            Assembler(**assemble_args).run()

        # Check for fragment length with single end
        assemble_args['single_end'] = True
        with self.assertRaisesRegexp(AssertionError, 'fragment length'):
            Assembler(**assemble_args).run()

    @staticmethod
    def generate_cell(cell_name):
        locus_names = ["TCRA", "TCRB"]
        imgt_seq_location = os.path.join(base_dir, "resources", "Mmus",
                                         "imgt_sequences")
        const_seq_file = os.path.join(base_dir, "resources", "Mmus",
                                      "constant_seqs.csv")
        out_dir = os.path.join(base_dir, "test_data", "results", "cell1")

        invariant_seqs = parse_invariant_cells(os.path.join(
            base_dir, "resources", "Mmus", "invariant_cells.json"))
        cell = parse_IgBLAST(
            'TCR', locus_names, out_dir, cell_name, imgt_seq_location, 'Mmus',
            'imgt', const_seq_file, invariant_seqs=invariant_seqs)
        return cell

    def test_parse_igblast(self):
        cell = self.generate_cell('cell1')
        assert not cell._check_is_empty(), "No Cell results"

    def test_align_salmon(self):
        cell = self.generate_cell('cell1')
        self.assembler.quantify(cell)

        # Look for recombinants in cell
        for receptor, locus_dict in six.iteritems(cell.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants:
                    for rec in recombinants:
                        print(rec.TPM)

    def xtest_invariant_seqs(self):
        cell = self.generate_cell('cell1')
        assert len(cell.invariant_seqs)
        all_vs = [seq['V'] for seq in cell.invariant_seqs]
        assert 'TRAV11' in all_vs
        assert not cell._check_if_inkt()


if __name__ == '__main__':
    unittest.main()
