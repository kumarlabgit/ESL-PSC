"""
Minimal FASTA reader used by the GUI.
"""

from typing import List, Tuple

def read_fasta(file_path: str) -> List[Tuple[str, str]]:
    """
    Parameters
    ----------
    file_path : str
        Path to the FASTA file.

    Returns
    -------
    List[Tuple[str, str]]
        List of (record_id, sequence) tuples in input order.
    """
    records: list[tuple[str, str]] = []
    with open(file_path, "r", encoding="utf-8") as handle:
        current_id: str | None = None
        current_seq: list[str] = []
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(current_seq)))
                current_id = line[1:].strip()
                current_seq.clear()
            else:
                current_seq.append(line)
        if current_id is not None:
            records.append((current_id, "".join(current_seq)))
    return records
