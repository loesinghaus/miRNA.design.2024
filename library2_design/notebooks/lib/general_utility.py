def intersection_of_lists(dictionary):
    # Convert the first list in the dictionary to a set
    intersection_set = set(dictionary[next(iter(dictionary))])
    
    # Iterate over the other lists in the dictionary and update the intersection set
    for key in dictionary:
        intersection_set.intersection_update(dictionary[key])

    return list(intersection_set)


