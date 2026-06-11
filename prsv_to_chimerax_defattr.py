#!/usr/bin/env python

import argparse
import csv
import os
import sys
from collections import defaultdict


ONE_LETTER_CODES = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert PRSV CSV values to a ChimeraX residue .defattr file."
    )
    parser.add_argument("prsv_csv", help="PRSV CSV with Res, PRSV, and Tubulin columns.")
    parser.add_argument("pdb_file", help="Target PDB model used to validate chain IDs and residues.")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output ChimeraX .defattr file. Default: PRSV CSV basename with .defattr extension.",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=1.0,
        help="Scale factor applied to every PRSV value. Default: %(default)s",
    )
    parser.add_argument(
        "--attribute-name",
        default="prsv",
        help="ChimeraX residue attribute name. Default: %(default)s",
    )
    parser.add_argument(
        "--model-number",
        default="1",
        help="ChimeraX model number used in residue selectors. Default: %(default)s",
    )
    parser.add_argument(
        "--alpha-chains",
        nargs="+",
        default=["A1", "A2"],
        help="Target chain IDs for Alpha rows. Default: A1 A2",
    )
    parser.add_argument(
        "--beta-chains",
        nargs="+",
        default=["B1", "B2"],
        help="Target chain IDs for Beta rows. Default: B1 B2",
    )
    parser.add_argument(
        "--include-missing",
        action="store_true",
        help="Write attributes even when a residue is not found in the target PDB.",
    )
    return parser.parse_args()


def output_path(prsv_csv, requested_output):
    if requested_output is not None:
        return requested_output

    base, _ = os.path.splitext(prsv_csv)
    return base + ".defattr"


def parse_pdb_residue(line):
    if not line.startswith("ATOM"):
        return None

    fields = line.split()
    if len(fields) < 5:
        return None

    residue_field = fields[3]

    if len(residue_field) > 3 and residue_field[:3] in ONE_LETTER_CODES:
        chain = residue_field[3:]
        residue_number = fields[4]
    elif len(fields) >= 6:
        chain = fields[4]
        residue_number = fields[5]
    else:
        return None

    try:
        residue_number = int(residue_number)
    except ValueError:
        return None

    return chain, residue_number


def pdb_residues_by_chain(pdb_file):
    residues = defaultdict(set)

    with open(pdb_file, "r") as handle:
        for line in handle:
            parsed = parse_pdb_residue(line)
            if parsed is None:
                continue

            chain, residue_number = parsed
            residues[chain].add(residue_number)

    return residues


def require_columns(reader, required_columns):
    missing = [column for column in required_columns if column not in reader.fieldnames]
    if missing:
        raise ValueError("PRSV CSV is missing required columns: %s" % ", ".join(missing))


def tubulin_chains(tubulin, alpha_chains, beta_chains):
    normalized = tubulin.strip().lower()

    if normalized == "alpha":
        return alpha_chains
    if normalized == "beta":
        return beta_chains

    raise ValueError("Unexpected Tubulin value %r. Expected Alpha or Beta." % tubulin)


def write_defattr(args):
    target_residues = pdb_residues_by_chain(args.pdb_file)
    requested_chains = set(args.alpha_chains + args.beta_chains)
    missing_chains = sorted(chain for chain in requested_chains if chain not in target_residues)

    if missing_chains:
        raise ValueError(
            "Target PDB is missing requested chain IDs: %s\nAvailable chain IDs: %s"
            % (", ".join(missing_chains), ", ".join(sorted(target_residues)))
        )

    output = output_path(args.prsv_csv, args.output)
    missing_residue_count = 0
    written_count = 0

    with open(args.prsv_csv, "r", newline="") as csv_handle, open(output, "w") as out_handle:
        reader = csv.DictReader(csv_handle)
        require_columns(reader, ["Res", "PRSV", "Tubulin"])

        out_handle.write("attribute: %s\n" % args.attribute_name)
        out_handle.write("recipient: residues\n")

        for row in reader:
            try:
                residue_number = int(float(row["Res"]))
                scaled_value = float(row["PRSV"]) * args.scale
            except ValueError as error:
                raise ValueError("Could not parse row as numeric Res/PRSV: %s" % row) from error

            chains = tubulin_chains(row["Tubulin"], args.alpha_chains, args.beta_chains)

            for chain in chains:
                residue_exists = residue_number in target_residues[chain]
                if not residue_exists:
                    missing_residue_count += 1
                    if not args.include_missing:
                        continue

                out_handle.write(
                    "\t#%s/%s:%s\t%s\n"
                    % (args.model_number, chain, residue_number, scaled_value)
                )
                written_count += 1

    if missing_residue_count:
        print(
            "Warning: %s chain/residue assignments were not found in the target PDB."
            % missing_residue_count,
            file=sys.stderr,
        )

    print("Wrote %s residue attributes to %s" % (written_count, output))


def main():
    args = parse_args()

    if args.scale <= 0:
        raise ValueError("--scale must be greater than 0")

    write_defattr(args)


if __name__ == "__main__":
    main()
