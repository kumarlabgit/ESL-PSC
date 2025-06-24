"""
Dialogs for displaying ESL-PSC analysis results.
"""
from __future__ import annotations
import os

# Use shared FASTA reader
from gui.core.fasta_io import read_fasta
import pandas as pd
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor, QBrush
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QSizePolicy, QTableWidget, QTableWidgetItem,
    QAbstractItemView, QTreeWidget, QTreeWidgetItem, QMessageBox,
    QLabel, QHeaderView
)
from PyQt6.QtSvgWidgets import QSvgWidget

from gui.ui.widgets.protein_map import ProteinMapWidget
from gui.ui.widgets.site_viewer import SiteViewer

# Keep references to open dialogs to prevent garbage collection
_open_dialogs = []


def _launch_site_viewer(gene: str, config, sites_path: str | None, parent=None) -> None:
    """Open the SiteViewer window for the given gene."""
    align_dir = getattr(config, "alignments_dir", "")
    if not align_dir:
        raise RuntimeError("Alignment directory not specified")
    align_path = os.path.join(align_dir, f"{gene}.fas")
    if not os.path.exists(align_path):
        raise FileNotFoundError(f"Alignment file not found: {align_path}")

    # Load alignment records using shared FASTA reader for reliability/consistency
    records = read_fasta(align_path)

    # Determine species groups from first response matrix
    response_dir = config.response_dir
    if not response_dir:
        base = os.path.basename(config.species_groups_file).replace(".txt", "")
        response_dir = os.path.join(config.output_dir, f"{base}_response_matrices")
    if not os.path.isdir(response_dir):
        raise FileNotFoundError(f"Response directory not found: {response_dir}")
    files = sorted([f for f in os.listdir(response_dir) if f.endswith('.txt')])
    if not files:
        raise FileNotFoundError("No response matrices found")
    first_matrix = os.path.join(response_dir, files[0])

    conv = []
    ctrl = []
    with open(first_matrix) as f:
        for line in f:
            sp, val = line.strip().split()[:2]
            if val == '1':
                conv.append(sp)
            else:
                ctrl.append(sp)

    # ─── Phenotype data ───────────────────────────────────────────────
    species_pheno_map: dict[str, int] | None = None
    pheno_name_map: dict[int, str] | None = None
    pheno_path = getattr(config, "species_phenotypes_file", "")
    if pheno_path and os.path.exists(pheno_path):
        try:
            species_pheno_map = {}
            with open(pheno_path) as fp:
                for line in fp:
                    parts = [p.strip() for p in line.strip().split(',') if p.strip()]
                    if len(parts) >= 2:
                        sp, val = parts[0], parts[1]
                        try:
                            species_pheno_map[sp] = int(val)
                        except ValueError:
                            # Skip malformed phenotype value
                            continue
        except Exception:
            species_pheno_map = None

    # Phenotype display names from config
    if hasattr(config, "pheno_name1") and hasattr(config, "pheno_name2"):
        pheno_name_map = {1: str(config.pheno_name1), -1: str(config.pheno_name2)}

    all_sites_info = None
    pss_map = {}
    if sites_path and os.path.exists(sites_path):
        try:
            df = pd.read_csv(sites_path)
            df = df[df["gene_name"] == gene]
            for _, row in df.iterrows():
                pos = int(row["position"]) - 1
                pss_map[pos] = float(row.get("pss", row.get("score", 0)))
        except Exception:
            pass

    # Launch as independent top-level window (no parent) so it opens separately
    viewer = SiteViewer(
        records,
        conv,
        ctrl,
        [],  # outgroup species (not derived here)
        gene,
        all_sites_info,
        False,
        pss_map,
        parent=None,
        species_pheno_map=species_pheno_map,
        pheno_name_map=pheno_name_map,
    )
    viewer.show()
    _open_dialogs.append(viewer)

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
        self.resize(680, 600)  

    @staticmethod
    def show_dialog(svg_path, parent=None):
        """Create, show, and store a reference to the dialog."""
        dialog = SpsPlotDialog(svg_path, parent)
        dialog.show()
        _open_dialogs.append(dialog)


class GeneRanksDialog(QDialog):
    """Dialog to display top gene ranks and optional selected sites."""

    def __init__(self, dataframe, config, sites_path=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Top Gene Ranks")
        layout = QVBoxLayout(self)

        self.config = config
        self.sites_path = sites_path

        df = dataframe

        # Clean and rename columns for readability
        rename_map = {
            'gene_name': 'Gene',
            'highest_gss': 'Highest GSS',
            'highest_ever_gss': 'Highest GSS',
            'best_rank': 'Best Rank',
            'best_ever_rank': 'Best Ever Rank',
            'num_selected_sites': 'Number of Selected Sites',
            'num_combos_ranked': 'Number of Combinations Ranked',
            'num_combos_ranked_top': 'Number Ranked in Top',
        }
        df = df.rename(columns={k: v for k, v in rename_map.items() if k in df.columns})

        # Round and format columns
        if 'Highest GSS' in df.columns:
            df['Highest GSS'] = df['Highest GSS'].round(5)
        # Ensure rank columns are numeric while safely handling placeholders like 'None'
        for col in ['Best Rank', 'Best Ever Rank']:
            if col in df.columns:
                # Convert values to numeric, setting invalid parsing (e.g. 'None', '') to NaN
                numeric_vals = pd.to_numeric(df[col], errors='coerce')
                # Replace NaN with empty string to keep the table display tidy
                df[col] = numeric_vals.apply(lambda x: '' if pd.isna(x) else int(x))

        self.has_sites = bool(sites_path and os.path.exists(sites_path))

        msg = "Double-click a row to open the Site Viewer."
        if self.has_sites:
            msg += " Use the expand/collapse arrows to reveal selected sites for each gene."
        help_label = QLabel(msg)
        help_label.setWordWrap(True)
        layout.addWidget(help_label)

        if self.has_sites:
            self._init_tree_view(layout, df)
        else:
            self._init_table_view(layout, df)

        # Set a wider initial size to better accommodate the content
        self.resize(1200, 800)
        
        # Ensure the window is scrollable if content is too wide
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        
        # Set a minimum size to prevent the window from being too small
        self.setMinimumSize(1000, 600)

    def _open_site_viewer(self, row_or_item, _column: int) -> None:
        """Launch :class:`SiteViewer` for the selected gene."""
        # When using the tree view the signal provides a ``QTreeWidgetItem``
        # whereas the table view provides the row index.  Handle both cases.
        if self.has_sites:
            # If called from ``QTreeWidget.itemDoubleClicked`` we are passed the
            # item directly.  The older implementation expected a row index,
            # which caused a crash because ``topLevelItem`` requires an int.  To
            # remain backwards compatible we accept either form.
            if isinstance(row_or_item, QTreeWidgetItem):
                item = row_or_item
            else:
                item = self.tree.topLevelItem(int(row_or_item)) if row_or_item is not None else None
            # Only act on top-level gene rows, not the child site rows.
            if item is None or item.parent() is not None:
                return
            gene = item.text(0)
        else:
            row = int(row_or_item)
            gene_item = self.table.item(row, self.gene_col)
            if gene_item is None:
                return
            gene = gene_item.text()
        try:
            _launch_site_viewer(gene, self.config, self.sites_path, parent=self)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to open Site Viewer:\n{e}")

    def _init_table_view(self, layout: QVBoxLayout, df: pd.DataFrame) -> None:
        """Initialize the simple table view (no selected sites)."""
        self.table = QTableWidget(len(df.index), len(df.columns))
        self.table.setHorizontalHeaderLabels(list(df.columns))
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSortingEnabled(False)
        self.gene_col = list(df.columns).index('Gene') if 'Gene' in df.columns else 0

        for r_idx, (_, row) in enumerate(df.iterrows()):
            for c_idx, value in enumerate(row):
                item = QTableWidgetItem(str(value))
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r_idx, c_idx, item)

        self.table.resizeColumnsToContents()
        layout.addWidget(self.table)
        self.table.cellDoubleClicked.connect(self._open_site_viewer)

    def _init_tree_view(self, layout: QVBoxLayout, df: pd.DataFrame) -> None:
        """Initialize the tree view with selected sites information."""
        try:
            sites_df = pd.read_csv(self.sites_path)
        except Exception:
            # Fallback to table view if parsing fails
            self.has_sites = False
            self._init_table_view(layout, df)
            return

        self.tree = QTreeWidget()

        rank_columns = [c for c in df.columns if c != 'Gene']
        headers = ['Gene', 'Length', 'Map', 'Position', 'Max Position Sparsity Score'] + rank_columns
        self.tree.setColumnCount(len(headers))
        self.tree.setHeaderLabels(headers)
        self.tree.setUniformRowHeights(True)
        self.tree.header().setSectionResizeMode(2, QHeaderView.ResizeMode.Fixed)
        self.tree.setColumnWidth(2, 520)

        align_dir = getattr(self.config, 'alignments_dir', '')

        for _, row in df.iterrows():
            gene = row['Gene']
            gene_sites = sites_df[sites_df['gene_name'] == gene]
            gene_sites = gene_sites[gene_sites['pss'] != 0] if 'pss' in gene_sites.columns else gene_sites
            positions = gene_sites['position'].tolist()
            scores = (
                gene_sites['pss'].tolist()
                if 'pss' in gene_sites.columns and not gene_sites['pss'].isnull().all()
                else gene_sites['score'].tolist() if 'score' in gene_sites.columns else []
            )
            # Pre-compute score range for consistent coloring
            min_score = min(scores) if scores else None
            max_score = max(scores) if scores else None
            if min_score is not None and max_score == min_score:
                max_score += 1e-9
            length = _get_alignment_length(gene, align_dir) or (max(positions) if positions else 1)

            gene_values = [gene, str(length), '', '', ''] + [str(row[col]) for col in rank_columns]
            parent_item = QTreeWidgetItem(gene_values)
            self.tree.addTopLevelItem(parent_item)
            self.tree.setItemWidget(parent_item, 2, ProteinMapWidget(length, positions, scores))
            parent_item.setExpanded(False)

            for _, srow in gene_sites.iterrows():
                pos = str(srow['position'])
                pss_val = srow['pss'] if 'pss' in srow else srow.get('score', 0)
                if pss_val == 0:
                    continue
                pss = f"{pss_val:.5f}"
                child_vals = ['', '', '', pos, pss] + ['' for _ in rank_columns]
                child = QTreeWidgetItem(child_vals)
                # Apply the same color scaling used by ProteinMapWidget
                if min_score is not None and max_score is not None:
                    rel = (pss_val - min_score) / (max_score - min_score)
                    brightness = 0.2 + 0.8 * rel
                else:
                    brightness = 1.0
                base_color = QColor("#4CAF50")
                r = int(base_color.red() * brightness)
                g = int(base_color.green() * brightness)
                b = int(base_color.blue() * brightness)
                brush = QBrush(QColor(r, g, b))
                child.setForeground(4, brush)
                parent_item.addChild(child)

        self.tree.resizeColumnToContents(0)
        self.tree.resizeColumnToContents(1)
        self.tree.resizeColumnToContents(3)
        self.tree.resizeColumnToContents(4)

        layout.addWidget(self.tree)
        self.tree.itemDoubleClicked.connect(self._open_site_viewer)

    @staticmethod
    def show_dialog(csv_path, config, sites_path=None, parent=None):
        """Load data, then create, show, and store a reference to the dialog."""
        try:
            # Read file as-is and keep the existing row order (already sorted by the CLI).
            df = pd.read_csv(csv_path).head(200)
        except Exception as e:
            QMessageBox.critical(parent, "Error", f"Failed to read gene ranks file:\n{csv_path}\n\n{e}")
            return

        dialog = GeneRanksDialog(df, config, sites_path, parent)
        dialog.show()
        _open_dialogs.append(dialog)


class SelectedSitesDialog(QDialog):
    """(Deprecated) Dialog to display selected sites grouped by gene."""
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
            scores = (
                gene_data["pss"].tolist()
                if "pss" in gene_data.columns and not gene_data["pss"].isnull().all()
                else gene_data["score"].tolist() if "score" in gene_data.columns else []
            )
            length = _get_alignment_length(gene, alignments_dir) or (max(positions) if positions else 1)

            parent_item = QTreeWidgetItem([gene, str(length), "", "", ""])
            tree.addTopLevelItem(parent_item)
            tree.setItemWidget(parent_item, 2, ProteinMapWidget(length, positions, scores))
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