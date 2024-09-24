import numpy as np
#import logomaker
import pandas as pd
import matplotlib.pyplot as plt
from typing import Union

# alphabets

def get_alphabet(alph: str = "DNA"):
    DNA_alphabet = ["A", "C", "G", "T"]
    RNA_alphabet = ["A", "C", "G", "U"]

    if alph == "DNA":
        return DNA_alphabet
    elif alph == "RNA":
        return RNA_alphabet
    else:
        raise ValueError("Invalid alphabet. Please choose 'DNA' or 'RNA'.")


def one_hot_encoding(sequences: Union[list, str], alph: str = "DNA"):
    """Expects a list of sequences or a single string.

    Returns a one-hot encoded numpy array."""
    alphabet = get_alphabet(alph)

    # create mapping
    char_to_int = {c: i for i, c in enumerate(alphabet)}
    unpack_flag = False
    if not isinstance(sequences, list):
        unpack_flag = True
        sequences = [sequences]
    one_hot = np.zeros((len(sequences), len(sequences[0]), len(alphabet)))
    for index, sequence in enumerate(sequences):
        # convert to integer encoding
        integer_encoded = [char_to_int[c] for c in sequence]
        # one hot encoding
        one_hot[index, :, :] = np.eye(len(alphabet))[integer_encoded]
    if unpack_flag:
        one_hot = one_hot.squeeze(axis=0)
    return one_hot


def one_hot_decoding(one_hot: np.array, alph: str = "DNA"):
    """ Expects a numpy array in the form [sample, sequence, letter].

    Returns a list of sequences."""
    alphabet = get_alphabet(alph)

    int_array = np.argmax(one_hot, axis=-1)
    # int_to_char = {i: c for i, c in enumerate(alphabet)}
    # keys = np.array(list(int_to_char.keys()))
    # values = np.array(list(int_to_char.values()))
    # mapping_array = np.zeros(keys.max()+1, dtype=values.dtype)
    # mapping_array[keys] = values
    int_to_char = np.array(alphabet)
    char_array = int_to_char[int_array]
    string_list = [''.join(row) for row in char_array]
    return string_list


def calc_PPM(sequences: list, alph: str = "DNA", pseudocounts: float = 0.01):
    alphabet = get_alphabet(alph)

    # count calculation
    one_hot = one_hot_encoding(sequences, alph)
    counts = np.sum(one_hot, axis=0)
    # add pseudocounts
    counts += pseudocounts
    # frequences
    freqs = counts/one_hot.shape[0]
    return freqs


def calc_PWM(PPM: np.array, alph: str = "DNA"):
    return np.log2(PPM/0.25)


def calc_information_content(PPM: np.array, alph: str = "DNA"):
    alphabet = get_alphabet(alph)
    IC_total = np.log2(len(alphabet))
    uncertainty_per_position = -np.sum(PPM*np.log2(PPM), axis=1, keepdims=True)
    IC_final = IC_total - uncertainty_per_position
    info_per_letter = PPM*IC_final

    return info_per_letter


# def plot_logo(matrix: np.array, alph: str = "DNA"):
#     alphabet = get_alphabet(alph)
#     logo_data = pd.DataFrame(matrix, columns=alphabet)
#     logo = logomaker.Logo(logo_data,
#                           shade_below=.5,
#                           fade_below=.5,
#                           font_name='Arial Rounded MT Bold')
#     print(np.max(logo_data))
#     logo.ax.set_ylim([0, np.max(matrix)])
#     plt.show()

def reverse_complement(seq: str, alph: str = "RNA"):
    if alph == "DNA":
        complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    elif alph == "RNA":
        complements = {"A": "U", "C": "G", "G": "C", "U": "A"}
    else:
        raise ValueError("Invalid alphabet. Please choose 'DNA' or 'RNA'.")

    seq = seq[::-1]
    seq = "".join([complements[c] for c in seq])

    return seq

def complement(seq: str, alph: str = "RNA"):
    if alph == "DNA":
        complements = {"A": "T", "C": "G", "G": "C", "T": "A"}
    elif alph == "RNA":
        complements = {"A": "U", "C": "G", "G": "C", "U": "A"}
    else:
        raise ValueError("Invalid alphabet. Please choose 'DNA' or 'RNA'.")

    seq = "".join([complements[c] for c in seq])

    return seq

def one_point_mutation(primer):
    nucleotides = ["A", "T", "G", "C"]
    for index, base in enumerate(primer):
        for nt in nucleotides:
            if base != nt:
                yield primer[:index] + nt + primer[index+1:]