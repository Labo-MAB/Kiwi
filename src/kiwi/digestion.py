"""
This script define the class Digestion which perform the in silico digestion and
check the uniqueness of peptide sequences.

It contains also the lines for command line usage.
"""

import re
import numpy as np

from kiwi.reader import read_fasta
from kiwi.enzyme import cleavage_site
from kiwi.mass import amino_acid, atom


class Experiment:
    """Instantiate this class to perform in silico digestion and find unique
    peptide"""

    def __init__(
        self,
        fasta,
        enzyme="trypsin",
        min_length=7,
        max_length=None,
        max_miscleavages=1,
        max_mass=4600,
    ):

        if isinstance(fasta, str):
            self.fasta = fasta
            self.proteins = read_fasta(fasta)
        elif isinstance(fasta, dict):
            self.fasta = ""
            self.proteins = fasta
        else:
            raise TypeError(
                f"Expected ``fasta`` of type str or dict. Received: {type(fasta)}"
            )
        self.params = {
            "enzyme": enzyme,
            "min_length": min_length,
            "max_length": max_length,
            "max_miscleavages": max_miscleavages,
            "max_mass": max_mass,
            "clip_nterm_m": True,
        }
        self.peptides = {id: [] for id in self.proteins}
        self.outdir = f"{fasta}_digestedPeptides.csv"
        self.is_unique = {}

    def set_enzyme(self, enzyme_name):
        self.params["enzyme"] = enzyme_name.lower()

    def set_min_length(self, minimum):
        self.params["min_length"] = int(minimum)

    def set_max_length(self, maximum):
        self.params["max_length"] = int(maximum)

    def set_max_miscleavages(self, miscleavages):
        self.params["max_miscleavages"] = int(miscleavages)

    def set_max_mass(self, _mass):
        self.params["max_mass"] = float(_mass)

    def set_clip_nterm_m(self, bool_):
        self.params["clip_nterm_m"] = bool_

    def set_outdir(self, dir_):
        self.outdir = dir_

    def _is_cleaved(self):
        return (
            len([pep for peptides in self.peptides.values() for pep in peptides]) != 0
        )

    def _get_mass(self, seq):
        """Sum the mass of all residues and add the mass of H (N-terminal) and OH
        (C-terminus) to obtain the sequence mass"""
        _mass = sum((amino_acid[aa] for aa in seq)) + atom["H"] * 2 + atom["O"]
        return round(_mass, 4)

    def _join_sequences(self, base_sequences, miss):
        """Iterate through base sequences (fully digested peptide) and join miss+1
        adjacent sequences"""
        num_sequences = len(base_sequences) - miss  # Number of sequences to loop on
        return ["".join(base_sequences[n : n + miss + 1]) for n in range(num_sequences)]

    def _is_good_peptide(self, seq):
        """Filter sequences with given threshold(s)"""
        if self.params["min_length"] is not None and self.params["min_length"] > len(
            seq
        ):
            return False
        if self.params["max_length"] is not None and self.params["max_length"] < len(
            seq
        ):
            return False
        try:
            if self.params["max_mass"] is not None and self.params[
                "max_mass"
            ] < self._get_mass(seq):
                return False
        # Catch non-amino acid characters (such as X)
        except KeyError:
            return False

        return True

    def _format_line(self, peptide, id_, include_unique_flag=False):
        if include_unique_flag:
            uniqueness = "unique" if self.is_unique[peptide] else "non-unique"
            return f"{peptide}, {id_}, {self._get_mass(peptide)}, {uniqueness}\n"

        return f"{peptide}, {id_}, {self._get_mass(peptide)}\n"

    def write(self, outdir=None):
        """Write results to a csv file"""
        if not self._is_cleaved():
            raise ValueError(
                "No sequences found. Please use cleave_proteins() to generate sequences"
            )

        if outdir is not None:
            self.set_outdir(outdir)

        include_unique_flag = len(self.is_unique) != 0

        with open(self.outdir, "w", encoding="UTF-8") as out_file:
            for id_ in self.peptides:
                for peptide in self.peptides[id_]:
                    out_file.write(self._format_line(peptide, id_, include_unique_flag))

    def cleave_proteins(self):
        """Perform an in silico digestion of proteins with previously selected
        or default parameters"""
        regexp = cleavage_site[self.params["enzyme"]]

        for id_, protein_sequences in self.proteins.items():
            # split protein into perfectly digested sequences
            base_sequences = re.split(regexp, protein_sequences)
            miss = 0
            while miss <= self.params["max_miscleavages"]:
                # Join base sequences in function of actual miscleavage value
                all_peptides = self._join_sequences(base_sequences, miss)

                if all_peptides == []:
                    miss += 1
                    continue

                # Also consider N-terminal peptide without its methionine
                if self.params["clip_nterm_m"]:
                    first_pep = all_peptides[0]
                    all_peptides.append(re.sub(r"^[Mm]", "", first_pep))

                # Filter all peptide with given threshold(s)
                if len(all_peptides) > 0:
                    self.peptides[id_].extend(
                        [pep for pep in all_peptides if self._is_good_peptide(pep)]
                    )

                miss += 1

    def check_peptide_uniqueness(self):
        """Check if peptide sequences are unique in the whole fasta file (i.e.
        has no identical match in other proteins)"""

        # Get all peptides. Use set() to remove peptide duplicates from same protein (if any)
        all_peptides = [
            peptide for peptides in self.peptides.values() for peptide in set(peptides)
        ]

        # Find frequency of each peptide (does not account peptide duplicates from same protein)
        peptides, count = np.unique(all_peptides, return_counts=True)
        peptides_frequency = dict(zip(peptides, count))

        self.is_unique = {
            pep: (peptides_frequency[pep] == 1)
            for peptides in self.peptides.values()
            for pep in peptides
        }
