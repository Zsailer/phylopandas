import pandas as pd

COLORS = {
    "A": "#999999",
    "C": "#FFFF00",
    "D": "#33CCFF",
    "E": "#33CCFF",
    "F": "#FF9900",
    "G": "#999999",
    "H": "#CC0000",
    "I": "#FF6666",
    "K": "#CC0000",
    "L": "#FF6666",
    "M": "#3366FF",
    "N": "#3366FF",
    "P": "#CC33CC",
    "Q": "#3366FF",
    "R": "#CC0000",
    "S": "#999999",
    "T": "#3366FF",
    "V": "#FF6666",
    "W": "#FF9900",
    "Y": "#FF9900",
    "-": "#F5F5F5"
}

def get_alignment(df, index_col='id', apply_colors=False, gap_character='-'):
    """Build an alignment DataFrame.
    
    Rows are indexed by ID. Columns are labeled by position number.
    """
    sequences = df.sequence

    # Get dimensions
    nseq = len(sequences)
    length = max([len(s) for s in sequences])

    # Build DataFrame
    columns = range(length)
    index = df[index_col]
    alignment = pd.DataFrame(index=index, columns=columns)
    alignment.fillna(value=gap_character)

    # Replace Alignment with residues.
    for i, idx in enumerate(index):
        alignment.loc[idx] = list(sequences.iloc[i])

    # Apply residue colors for display.
    if apply_colors:
        colormap = lambda char: "background-color: %s;" % COLORS[char]
        alignment = alignment.style.applymap(colormap)
    return alignment