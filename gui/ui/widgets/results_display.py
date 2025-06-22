"""
Dialogs for displaying ESL-PSC analysis results.
"""
import os
import pandas as pd
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QSizePolicy, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QTreeWidget, QTreeWidgetItem, QMessageBox
)
from PyQt6.QtSvgWidgets import QSvgWidget

from gui.ui.widgets.protein_map import ProteinMapWidget

# Keep references to open dialogs to prevent garbage collection
_open_dialogs = []

def _get_alignment_length(gene_name: str, alignments_dir: str):
    """Return the length of the alignment for the given gene."""
    if not alignments_dir:
        return None
    path = os.path.join(alignments_dir, f"{gene_name}.fas")
    try:
        with open(path) as handle:
            for line in handle:
                if not line.startswith(">"):
                    return len(line.strip())
    except OSError:
        return None
    return None


class SpsPlotDialog(QDialog):
    """Dialog to display the SVG SPS plot."""
    def __init__(self, svg_path, parent=None):
        super().__init__(parent)
        self.setWindowTitle("SPS Prediction Plot")
        layout = QVBoxLayout(self)
        svg_widget = QSvgWidget(svg_path)
        svg_widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(svg_widget)
        self.resize(800, 600)

    @staticmethod
    def show_dialog(svg_path, parent=None):
        """Create, show, and store a reference to the dialog."""
        dialog = SpsPlotDialog(svg_path, parent)
        dialog.show()
        _open_dialogs.append(dialog)


class GeneRanksDialog(QDialog):
    """Dialog to display top gene ranks from a pandas DataFrame."""
    def __init__(self, dataframe, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Top Gene Ranks")
        layout = QVBoxLayout(self)

        df = dataframe

        # Clean and rename columns
        rename_map = {
            'gene_name': 'Gene',
            'highest_gss': 'Highest GSS',
            'highest_ever_gss': 'Highest GSS',
            'best_rank': 'Best Rank',
            'best_ever_rank': 'Best Ever Rank',
            'num_selected_sites': 'Number of Selected Sites',
            'num_combos_ranked': 'Number of Combinations Ranked',
            'num_combos_ranked_top': 'Number Ranked in Top'
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        # Round GSS columns
        if 'Highest GSS' in df.columns:
            df['Highest GSS'] = df['Highest GSS'].round(5)
        # Ensure rank columns are ints
        for col in ['Best Rank', 'Best Ever Rank']:
            if col in df.columns:
                df[col] = df[col].fillna('').apply(lambda x: int(float(x)) if str(x).strip() != '' else '')

        table = QTableWidget(len(df.index), len(df.columns))
        table.setHorizontalHeaderLabels(list(df.columns))
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        # Disable sorting so initial display matches CSV order.
        table.setSortingEnabled(False)

        for r_idx, (_, row) in enumerate(df.iterrows()):
            for c_idx, value in enumerate(row):
                item = QTableWidgetItem(str(value))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                table.setItem(r_idx, c_idx, item)

        table.resizeColumnsToContents()
        # Leave sorting disabled by default; users can manually enable sorting by clicking headers if desired.
        layout.addWidget(table)
        self.resize(900, 600)

    @staticmethod
    def show_dialog(csv_path, parent=None):
        """Load data, then create, show, and store a reference to the dialog."""
        try:
            # Read file as-is and keep the existing row order (already sorted by the CLI).
            df = pd.read_csv(csv_path).head(200)
        except Exception as e:
            QMessageBox.critical(parent, "Error", f"Failed to read gene ranks file:\n{csv_path}\n\n{e}")
            return

        dialog = GeneRanksDialog(df, parent)
        dialog.show()
        _open_dialogs.append(dialog)


class SelectedSitesDialog(QDialog):
    """Dialog to display selected sites grouped by gene."""
    def __init__(self, dataframe, gene_order, alignments_dir, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Selected Sites")
        layout = QVBoxLayout(self)
        
        df = dataframe

        tree = QTreeWidget()
        tree.setColumnCount(5)
        tree.setHeaderLabels(["Gene", "Length", "Map", "Position", "Position Sparsity Score"])
        tree.setUniformRowHeights(True)

        for gene in gene_order:
            gene_data = df[df["gene_name"] == gene]
            # Filter out rows with zero PSS
            gene_data = gene_data[gene_data["pss"] != 0] if "pss" in gene_data.columns else gene_data
            if gene_data.empty:
                continue

            positions = gene_data["position"].tolist()
            length = _get_alignment_length(gene, alignments_dir) or (max(positions) if positions else 1)

            parent_item = QTreeWidgetItem([gene, str(length), "", "", ""])
            tree.addTopLevelItem(parent_item)
            tree.setItemWidget(parent_item, 2, ProteinMapWidget(length, positions))
            parent_item.setExpanded(False)

            for _, row in gene_data.iterrows():
                pos = str(row["position"])
                pss_val = row["pss"] if "pss" in row else row.get("score", 0)
                if pss_val == 0:
                    continue  # skip zero PSS
                pss = f"{pss_val:.5f}"
                child = QTreeWidgetItem(["", "", "", pos, pss])
                parent_item.addChild(child)
        
        tree.resizeColumnToContents(0)
        tree.resizeColumnToContents(1)
        tree.resizeColumnToContents(3)
        tree.resizeColumnToContents(4)

        layout.addWidget(tree)
        self.resize(800, 600)

    @staticmethod
    def show_dialog(csv_path, gene_ranks_path, alignments_dir, parent=None):
        """Load data, then create, show, and store a reference to the dialog."""
        try:
            df = pd.read_csv(csv_path)
        except Exception as e:
            QMessageBox.critical(parent, "Error", f"Failed to read selected sites file:\n{csv_path}\n\n{e}")
            return

        # Determine gene display order
        gene_order = []
        if gene_ranks_path and os.path.exists(gene_ranks_path):
            try:
                ranks_df = pd.read_csv(gene_ranks_path)
                gene_order = list(ranks_df["gene_name"].head(200))
            except Exception:
                pass  # Fallback to order in sites file

        if not gene_order:
            gene_order = list(dict.fromkeys(df["gene_name"]))[:200]

        dialog = SelectedSitesDialog(df, gene_order, alignments_dir, parent)
        dialog.show()
        _open_dialogs.append(dialog)