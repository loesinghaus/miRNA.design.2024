import matplotlib.cm as cm

tab10 = cm.get_cmap('tab10')
tab20b = cm.get_cmap('tab20b')

cell_line_colors = {
    "HEK293T": tab10(0),
    "HeLa": tab10(4),
    "MCF7": tab10(2),
    "A549": tab10(3),
    "HaCaT": tab10(1), 
    "HUH7": tab10(5),
    "PC3": tab10(6),
    "JEG3": tab10(7),
    "Tera1": tab10(8),
    "SKNSH": tab10(9), 
    "K562": tab20b(5) 
}

cell_line_symbols = {
    "HEK293T": "o",  # Circle
    "HeLa": "D",  # Diamond
    "MCF7": "v",  # Triangle down
    "A549": "^",  # Triangle up
    "HaCaT": "<",  # Triangle left
    "HUH7": ">",  # Triangle right
    "PC3": "p",  # Pentagon
    "JEG3": "*",  # Star
    "Tera1": "h",  # Hexagon
    "SKNSH": "s",  # Square
    "K562": "+"    # Plus
}