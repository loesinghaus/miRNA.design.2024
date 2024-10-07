import numpy as np
from typing import Union
from .NA_sequence_utilities import get_alphabet

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
    """Expects a numpy array in the form [sample, sequence, letter].

    Returns a list of NA sequences."""
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
    """Expects a list of sequences.
    Returns a position probability matrix (PPM)."""

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
    """Expects a position probability matrix (PPM).
    Returns a position weight matrix (PWM)."""
    return np.log2(PPM/0.25)

def calc_information_content(PPM: np.array, alph: str = "DNA"):
    """Expects a position probability matrix (PPM).
    Returns the information content of the sequence per letter."""
    alphabet = get_alphabet(alph)
    IC_total = np.log2(len(alphabet))
    uncertainty_per_position = -np.sum(PPM*np.log2(PPM), axis=1, keepdims=True)
    IC_final = IC_total - uncertainty_per_position
    info_per_letter = PPM*IC_final

    return info_per_letter