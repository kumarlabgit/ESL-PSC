# GUI Quick Start

This guide is for first-time users who want the desktop app and do not want to set up the command line.

## Download and install

Use the GUI download for your computer:

- Mac with Apple Silicon: [`ESL-PSC-v2.5.1-macOS.dmg`](https://github.com/John-Allard/ESL-PSC/releases/download/v2.5.1/ESL-PSC-v2.5.1-macOS.dmg)
- Windows: [`ESL-PSC-v2.5.1-Windows.zip`](https://github.com/John-Allard/ESL-PSC/releases/download/v2.5.1/ESL-PSC-v2.5.1-Windows.zip)
- Ubuntu or Debian Linux on Intel/AMD 64-bit computers: [`esl-psc-gui_2.5.1_amd64.deb`](https://github.com/John-Allard/ESL-PSC/releases/download/v2.5.1/esl-psc-gui_2.5.1_amd64.deb)


Install the app:

### Mac

1. Download the macOS `.dmg` file. The current macOS GUI release is built for Apple Silicon Macs.
2. Open it.
3. Drag `ESL-PSC.app` into your `Applications` folder.
4. Open `ESL-PSC` from `Applications`, Launchpad, or Spotlight.

### Windows

1. Download the `.zip` file.
2. Open it and choose **Extract All**.
3. Open the extracted folder.
4. Double-click `ESL-PSC.exe`.

Windows may warn that the app is unsigned. If that happens, click **More info**, then **Run anyway**.

### Ubuntu or Debian Linux

1. Download the `.deb` file.
2. Double-click it and install it with your normal software installer.
3. Open `ESL-PSC` from your app menu.

The `.deb` file is for 64-bit Ubuntu, Debian, and closely related Debian-based Linux distributions. Other Linux distributions or ARM Linux systems should use the source-install instructions in the main README.

If you prefer the terminal on Linux, you can also install it with:

```bash
sudo dpkg -i esl-psc-gui_2.5.1_amd64.deb
```

## Getting started

Before you open the program, gather these files:

- a folder of alignment files in FASTA format
- a Newick tree file
- optionally, a CSV file with phenotype annotations

The phenotype CSV is helpful for tree coloring and for prediction summaries, but it is not required to start.

## Input file formats

### Alignment files

The alignment input should be a folder containing one FASTA file per gene or protein. Each file should contain aligned sequences for homologous species, and the species names in the FASTA headers should match the species names used in the tree and phenotype files.

ESL-PSC expects each sequence to be written on one line after its header. For example:

```text
>Homo_sapiens
MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGLGWAGWLLSPRGSRPSWGPTDPRRRSRHWV
>Mus_musculus
MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGLGWAGWLLSPRGSRPSWGPTDPRRRSRHWV
```

If your FASTA files wrap sequences across multiple lines, the GUI can convert them to the required two-line FASTA format before running the analysis.

### Newick tree file

The tree file should be a rooted Newick tree containing the species you want to use for contrast-pair selection. Branch lengths are allowed. Species labels should match the alignment headers and phenotype names.

Example:

```text
((Homo_sapiens:0.1,Pan_troglodytes:0.1):0.2,Mus_musculus:0.3);
```

### Phenotype file

The optional phenotype file should be a CSV-style text file with one species per line. For binary traits, use `1` for the convergent or focal trait state and `-1` for the control state:

```text
Homo_sapiens,1
Pan_troglodytes,1
Mus_musculus,-1
```

For continuous traits, the second column can be a numeric value:

```text
Homo_sapiens,1.82
Pan_troglodytes,1.55
Mus_musculus,0.03
```

### Species groups file

The species groups file defines the paired-species contrasts used for model training. The easiest way to make this file is with the Tree Viewer, but it is also a plain text file. Each contrast pair is represented by two consecutive lines: the focal or convergent species group first, then the paired control species group.

```text
convergent_species_A,convergent_species_B
control_species_A,control_species_B
```

Multiple contrast pairs are written as additional two-line blocks. Species names should match the alignment headers.

### 1. Fill in the Input page

![Annotated ESL-PSC input page](../images/quickstart-input.png)

1. Choose the folder that contains your alignments.
2. Click **Create a Species Groups File Using a Newick Tree**. A file chooser will open so you can pick your Newick tree.
3. After you save a species groups file from the tree viewer, it will appear here.
4. If you have a phenotype CSV file, choose it here. This is optional.
5. Click **Next**.

### 2. Build a species groups file from the tree

![Annotated tree viewer](../images/quickstart-tree.png)

1. If you have a phenotype CSV file, load it here so the tree can color the species labels.
2. Click **Auto Select Contrast Pairs** to let the program choose a set of valid pairs from the tree.
3. Click **Save Species Groups** and choose where to save the file.

After you save the species groups file, close the tree viewer and return to the main wizard. The species groups path should now be filled in on the Input page automatically.

### 3. Use the defaults on the Parameters page

![Annotated parameters page](../images/quickstart-parameters.png)

1. Choose an output folder for the results.
2. Choose which outputs to generate. **Gene ranks only** skips the species-predictions file and is much faster for proteome-scale analyses whose main goal is to identify candidate convergent genes rather than predict trait values for species.
3. For a first run, leave the rest of this page at the defaults.
4. Click **Next**.

### 4. Skip the Command page

The next page shows the exact terminal command that the GUI is going to run.

Most users do not need to do anything here. You can simply click **Next** again.

### 5. Run the analysis

On the final page:

1. Click **Run Analysis**.
2. Watch the progress in the built-in terminal panel.
3. When the run finishes, move on to the results page to inspect the ranked genes, predictions, and selected sites.

## First practice run

If you installed the desktop GUI from a release download, download and unzip the
separate [`ESL-PSC-v2.5.1-demo-data.zip`](https://github.com/John-Allard/ESL-PSC/releases/download/v2.5.1/ESL-PSC-v2.5.1-demo-data.zip)
archive first.

For the photosynthesis practice run, use:

- alignments folder: `ESL-PSC-demo-data/photosynthesis/alignments`
- tree file: `ESL-PSC-demo-data/photosynthesis/photo_tree.nwk`
- phenotype file: `ESL-PSC-demo-data/photosynthesis/photo_species_phenotypes.txt`
