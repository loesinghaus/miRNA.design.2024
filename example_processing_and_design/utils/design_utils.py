from utils.transfer_functions import transfer_function_linear
import random
import pandas as pd
import numpy as np
import ast


def shuffle_index_tuple_order(df):
    """Takes a DataFrame with a MultiIndex of tuples and randomly shuffles the order of elements within each tuple index,
    while preserving the row order. This is done to make the crossover events more random.
    
    Args:
        df: DataFrame with MultiIndex of tuples
        
    Returns:
        DataFrame with shuffled tuple indices"""
    
    # Get the original index tuples
    index_tuples = list(df.index)
    
    # Shuffle each tuple while maintaining row order
    shuffled_tuples = []
    for tup in index_tuples:
        # Convert tuple to list for shuffling
        tup_list = list(tup)
        random.shuffle(tup_list)
        # Convert back to tuple
        shuffled_tuples.append(tuple(tup_list))
        
    # Create new MultiIndex with shuffled tuples
    new_index = pd.MultiIndex.from_tuples(shuffled_tuples, names=df.index.names)
    
    # Return DataFrame with new index
    return df.set_index(new_index)

def add_mirna_combs(mirna_expr, combs):
    """
    Takes a dataframe with microRNA expressions and tuples of combinations.
    Returns a dataframe with the added expression of the microRNAs in the constructs.
    
    Parameters:
    mirna_expr: dataframe with microRNA expression (index: microRNA names, columns: cell lines)
    combs: list of tuples (each tuple is a combination of miRNAs)
    
    Returns:
    DataFrame with combined expressions for each combination
    """
    # Create MultiIndex
    multiindex = pd.MultiIndex.from_tuples(combs, names=[f'miRNA{i+1}' for i in range(len(combs[0]))])
    
    # Pre-compute all sums
    sums = [mirna_expr.loc[list(comb),:].sum(axis=0).values for comb in combs]
    
    # Create DataFrame directly with the computed values
    added_expression = pd.DataFrame(
        sums,
        index=multiindex,
        columns=mirna_expr.columns
    )
    
    return added_expression

def calculate_mse(df, mse_target, loss_emphasis):
    """This function calculates the mean squared error of a design for a given target stability distribution.
    
    df: stability values for different cell lines in the columns and microRNAs in the rows.
    mse_target: dictionary with cell lines as keys and target stability as values.
    loss_emphasis: dictionary with cell lines as keys and loss emphasis as values."""
    
    mse = (df - mse_target)**2
    mse = mse.mul(loss_emphasis, axis=1)
    mse = mse.mean(axis=1)

    return mse

def calculate_linear_loss(df, target_vals, loss_emphasis):
    """This function calculates the linear loss of a design for a given target distribution.
    
    Parameters:
    df: DataFrame with stability values for different cell lines in columns and microRNAs in rows
    target_vals: dictionary with cell lines as keys and target stability as values
    loss_emphasis: dictionary with cell lines as keys and loss emphasis as values
    
    Returns:
    Series with linear loss values for each microRNA combination"""
    
    loss = abs(df - target_vals)
    loss = loss.mul(loss_emphasis, axis=1)
    loss = loss.mean(axis=1)

    return loss

def tsi(x):
    """This function calculates the TSI for a given combination of cell lines.
    The input is a numpy array with the expression of the cell lines in the columns
    and the microRNAs in the rows."""
    # check if the orientation of x is correct
    if x.shape[1] > 50:
        raise ValueError("The number of columns is very large. The columns are supposed to be cell lines. Is this the case?")

    # if x is not normalized yet, normalize it
    # THIS IS IMPORTANT
    x = x/x.max(axis=1, keepdims=True)

    tsi = np.sum(1-x, axis=1)/(x.shape[1]-1)
    
    return tsi

def calculate_ratio_fitness(df, mse_target, clip = 0):
    """Calculate the ratio of expression between high and low target cell lines.
    
    This function calculates the ratio between the mean expression in cell lines with high target values
    versus those with low target values. It can optionally clip low expression values in high target cell lines.
    
    Args:
        df (pd.DataFrame): DataFrame with expression values for different cell lines in columns and microRNAs in rows
        mse_target (dict): Dictionary mapping cell line names to their target expression values
        clip (float, optional): Threshold below which expression values in high target cell lines are set to 0. Defaults to 0.
    
    Returns:
        pd.Series: Ratio of mean expression in high target cell lines to mean expression in low target cell lines
        
    Raises:
        AssertionError: If there are more than two distinct target values in mse_target
    """
    
    # get the high and low values in the target
    high_val = max(mse_target.values())
    low_val = min(mse_target.values())

    # get the cell lines associated with the values
    low_cell_lines = [cell_line for cell_line in mse_target.keys() if mse_target[cell_line] == low_val]
    high_cell_lines = [cell_line for cell_line in mse_target.keys() if mse_target[cell_line] == high_val]

    # make sure that there are only two value types for ratio fitness
    assert len(low_cell_lines) + len(high_cell_lines) == len(mse_target)

    # we compare the mean expression of the low cell lines to the mean expression of the high cell lines
    expression_low = df[low_cell_lines].median(axis=1)
    expression_high = df[high_cell_lines].median(axis=1)

    # if the expression of the high cell lines is lower than clip, set it to 0
    # this essentially excludes designs that are not active in the high cell lines
    expression_high[expression_high < clip] = 0

    # calculate the ratio of the expression of the high cell lines to the expression of the low cell lines
    ratio = expression_high/expression_low

    return ratio

def calculate_penalties(pop, forbidden_combinations):
    """Calculate penalties for designs containing forbidden miRNA combinations.

    This function checks each design in a population for combinations of miRNAs that
    should not be used together, based on a dictionary of forbidden combinations.
    A large penalty is applied to any design containing forbidden combinations.

    Args:
        pop: List of tuples, where each tuple contains the miRNAs in a design.
             miRNAs can be strings representing miRNA names or "empty".
        forbidden_combinations: Dictionary mapping miRNA names to lists of other miRNAs
                             that should not be combined with that miRNA.

    Returns:
        List of numeric penalties, one per design in the population. A penalty of 1000
        is applied to designs with forbidden combinations, 0 to valid designs.

    Example:
        >>> pop = [('miR-223-3p', 'miR-150-5p', 'empty'), 
                  ('miR-146a-5p', 'miR-155-5p', 'empty')]
        >>> forbidden = {'miR-223-3p': ['miR-150-5p'], 
                        'miR-146a-5p': ['miR-155-5p']}
        >>> calculate_penalties(pop, forbidden)
        [1000, 1000]
    """
    penalty_number = 1000
    
    penalties = []
    for design in pop:
        has_penalty = False
        for design_mirna in design:
            other_mirnas = [mirna for mirna in design if mirna != design_mirna]
            if design_mirna == "empty":
                intersection = set()
            else:
                forbidden_mirnas_design = forbidden_combinations[design_mirna]
                intersection = set(forbidden_mirnas_design).intersection(set(other_mirnas))
                
            if len(intersection) > 0:
                # print(design_mirna, " ", intersection)
                penalties.append(penalty_number)
                has_penalty = True
                break
        if not has_penalty:
            penalties.append(0)

    return penalties

def calculate_fitness(population, expression, popt, loss_emphasis={}, mse_target=[], loss_type="mse",
    clip_min=0, clip_max=1, use_combination_penalty=True, forbidden_combinations=None):
    """This function calculates the fitness of a population of designs based on the projected stabiltiy and target.
    
    population: list of tuples of microRNAs
    expression: dataframe with microRNA expression (index: microRNA names, columns: cell lines)
    popt: tuple of parameters for the transfer function
    loss_emphasis: dictionary with cell lines as keys and loss emphasis as values
    mse_target: dictionary with cell lines as keys and target stabilities as values 
    loss_type: string with the type of loss function to use ("mse", "ratio", "linear")
    use_combination_penalty: boolean to use the combination penalty
    forbidden_combinations: dictionary with miRNA names as keys and lists of forbidden miRNAs as values
    clip_min: minimum expression value to clip
    clip_max: maximum expression value to clip

    Returns: Fitness values for the designs in the population.
    """
    if use_combination_penalty and forbidden_combinations is None:
        raise ValueError("Forbidden combinations are required for combination penalty.")
    
    # Calculate the stability levels for the designs in the population according to the additive model
    add_expr = add_mirna_combs(expression, population).apply(lambda x: transfer_function_linear(x, *popt))

    # if loss emphasis is empty, generate it
    if len(loss_emphasis) == 0:
        loss_emphasis = {cell_line: 1 for cell_line in add_expr.columns}

    if loss_type == "mse":
        fitness_val = 1./calculate_mse(add_expr, mse_target, loss_emphasis)
    elif loss_type == "linear":
        fitness_val = 1./calculate_linear_loss(add_expr, mse_target, loss_emphasis)
    elif loss_type == "ratio":
        fitness_val = calculate_ratio_fitness(add_expr, mse_target, clip=clip_min)
    else:
        raise(f"Loss type not implemented: {loss_type}")
    
    if use_combination_penalty:
        combination_penalty = calculate_penalties(population, forbidden_combinations)
    else:
        combination_penalty = 0
    
    fitness = fitness_val - combination_penalty

    return fitness

def evaluate_fitness(population, expression, popt, loss_emphasis={}, mse_target=[], loss_type="mse",
    clip_min=0, clip_max=1, use_combination_penalty=True, forbidden_combinations=None):
    """This function evaluates the fitness of a population of designs based on the projected stability and target.
    
    population: list of tuples of microRNAs
    expression: dataframe with microRNA expression (index: microRNA names, columns: cell lines)
    loss_emphasis: dictionary with cell lines as keys and loss emphasis as values
    mse_target: dictionary with cell lines as keys and target stabilities as values 
    loss_type: string with the type of loss function to use ("mse", "ratio", "linear")
    use_combination_penalty: boolean to use the combination penalty

    Returns: Fitness values for the designs in the population and the predicted stability values."""
    
    add_expr = add_mirna_combs(expression, population).apply(lambda x: transfer_function_linear(x, *popt))
    add_expr["quality"] = calculate_fitness(population, expression, popt, loss_emphasis, mse_target,
        loss_type, clip_min, clip_max, use_combination_penalty, forbidden_combinations)

    return add_expr

# -----------------------------------------------------------------------

def drop_duplicates(df):
    """ Drop all duplicate designs. Assumes that the indices are tuples of microRNAs. """
    sorted_idx = df.index.map(sorted)
    df['sorted_index'] = [tuple(i) for i in sorted_idx]
    duplicates = df.duplicated(subset='sorted_index', keep=False)
    df = df.drop_duplicates(subset='sorted_index', keep='first').drop(columns='sorted_index')
    return df, duplicates 

def select_parents(fitnesses, tournament_size=3):
    """Selects two parents from the population using tournament selection."""
    # Tournament selection
    parents = []

    for _ in range(2):  # Select two parents
        tournament = fitnesses.sample(tournament_size)
        # select the best individual
        winner = tournament.sort_values(ascending=False).index[0]
        parents.append(winner)

    return tuple(parents)

def crossover(parent1, parent2, n):
    """Performs single point crossover between two parent designs.
    
    Args:
        parent1: Tuple of microRNAs representing first parent design
        parent2: Tuple of microRNAs representing second parent design 
        n: Integer length of the designs
        
    Returns:
        Tuple containing the child design created by combining segments of both parents
    """
    # Single point crossover
    idx = random.randint(0, n-1)
    child = parent1[:idx] + parent2[idx:]
    return child

def mutate(child, miRNAs, n, mutation_rate=0.2):
    """Randomly mutates a child design by replacing one microRNA.
    
    Args:
        child: Tuple of microRNAs representing the design to mutate
        miRNAs: List of possible microRNAs to choose from
        n: Integer length of the design
        
    Returns:
        Tuple containing the mutated design. Has 20% chance of mutation,
        otherwise returns original design unchanged.
    """
    # Randomly replace one microRNA with another
    if random.random() < mutation_rate:  # 20% mutation rate
        idx = random.randint(0, n-1)
        new_mirna = random.choice(miRNAs)
        child = list(child)
        child[idx] = new_mirna
    return tuple(child)

# -----------------------------------------------------------------------

def determine_mirna_usage(df):
    usage_dict = {}
    used_mirnas = df.index.tolist()
    for design in used_mirnas:
        for mirna in design:
            if mirna in usage_dict:
                usage_dict[mirna] += 1
            else:
                usage_dict[mirna] = 1
    
    # sort dict by value
    usage_dict = {k: v for k, v in sorted(usage_dict.items(), key=lambda item: item[1], reverse=True)}
    return usage_dict

def count_mirnas_per_design(df):
    """Count the frequency of each miRNA across designs in a DataFrame.
    
    Args:
        df: DataFrame containing miRNA designs as index
        
    Returns:
        DataFrame containing counts of how many designs each miRNA appears in,
        sorted by frequency in descending order. If a miRNA appears multiple
        times in a single design, it is only counted once for that design.
    """
    combinations = df.index.tolist()

    mirna_count = {}
    for design in combinations:
        design_count = {}
        for mirna in design:
            if mirna in design_count:
                continue
            else:
                design_count[mirna] = 1
            if mirna in mirna_count:
                mirna_count[mirna] += 1
            else:
                mirna_count[mirna] = 1

    mirna_count_df = pd.DataFrame.from_dict(mirna_count, orient='index', columns=['count'])
    mirna_count_df = mirna_count_df.sort_values(by=['count'], ascending=False)            
    return mirna_count_df

def add_numbered_index(df, base_name):
    """Df is assumed to have a multiindex of microRNAs. First, convert the multi-index to columns.
    Then, add a column with the design number."""
    df = df.reset_index()
    df.index = [f"{base_name}_{i+1}" for i in range(len(df))]
    return df

# -----------------------------------------------------------------------

def generate_genetic_design(target, n_mirnas, mirnas, mirna_expression, popt, loss_emphasis={},
                            no_designs=10, generations=30, population_size=500, loss_type="mse", 
                            clip_min=0, clip_max=1, use_combination_penalty=True, forbidden_combinations=None):
    """Generate a single set of microRNA designs using a genetic algorithm approach.
    
    Args:
        target: Target expression values to optimize for
        n_mirnas: Number of microRNAs to include in each design
        mirnas: List of available microRNAs to choose from
        mirna_expression: DataFrame containing expression data for microRNAs
        popt: Tuple of parameters for the transfer function
        loss_emphasis: Dict specifying emphasis weights for loss calculation
        no_designs: Number of top designs to return
        generations: Number of generations to run the genetic algorithm
        population_size: Size of the population in each generation

        loss_type: Type of loss function to use ("mse" or other)
        use_combination_penalty: Whether to penalize forbbiden combinations of microRNAs
        
    Returns:
        DataFrame containing the top designs sorted by quality score, with columns for
        microRNA combinations and their performance metrics
    """
    
    # Initial population
    population = [tuple(random.choice(mirnas) for _ in range(n_mirnas)) for _ in range(population_size)]

    # Run the GA for a set number of generations
    for generation in range(generations):
        fitnesses = calculate_fitness(
            population=population,
            expression=mirna_expression,
            popt=popt,
            loss_emphasis=loss_emphasis,
            loss_type=loss_type,
            mse_target=target,
            clip_min=clip_min,
            clip_max=clip_max,
            use_combination_penalty=use_combination_penalty,
            forbidden_combinations=forbidden_combinations
            )
        new_population = []
        for _ in range(population_size):
            parent1, parent2 = select_parents(fitnesses)
            child = crossover(parent1, parent2, n_mirnas)
            child = mutate(child, mirnas, n_mirnas)
            new_population.append(child)
        population = new_population

    # Get the best designs
    designs = evaluate_fitness(
        population=population,
        expression=mirna_expression,
        loss_emphasis=loss_emphasis,
        mse_target=target,
        loss_type=loss_type,
        popt=popt,
        clip_min=clip_min,
        clip_max=clip_max,
        use_combination_penalty=use_combination_penalty,
        forbidden_combinations=forbidden_combinations
        )
    designs, _ = drop_duplicates(designs)
    designs = shuffle_index_tuple_order(designs)
    designs.sort_values(by=['quality'], ascending=False, inplace=True)
    designs = designs.head(no_designs)

    return designs

def generate_multiple_designs(target_values, designs_per_target, base_name, mirna_data, popt,
                    loss="mse", n_mirnas=5, loss_emphases={}, increase_diversity=1, clip_min=0, clip_max=1,
                    use_combination_penalty=True, forbidden_combinations=None):
    """Generate multiple microRNA designs targeting different MSE values.
    
    Args:
        target_values: List of target values to design for
        designs_per_target: Number of designs to generate per MSE target
        base_name: Base name for naming the designs
        mirna_data: DataFrame containing microRNA expression data
        popt: Tuple of parameters for the transfer function
        loss: Loss function to use ('mse' or other)
        n_mirnas: Number of microRNAs per design
        loss_emphases: Dict of emphasis values for different loss components
        increase_diversity: Factor to increase design diversity
        clip_min: Minimum value to clip expressions to
        clip_max: Maximum value to clip expressions to
        use_combination_penalty: Whether to penalize certain microRNA combinations
        forbidden_combinations: List of forbidden microRNA combinations
        
    Returns:
        DataFrame containing all generated designs with their properties
    """
    # List of microRNAs and their impacts  
    miRNAs = list(mirna_data.index)

    all_designs = []

    if len(loss_emphases) == 0:
        loss_emphases = [{} for i in range(len(target_values))]
        
    for i, mse_target in enumerate(target_values):
        print(f"Processing {base_name} {i+1}/{len(target_values)}")
        miRNAs_filter = miRNAs.copy()
        target_designs = []
        
        for _ in range(increase_diversity):
            designs = generate_genetic_design(
                target=mse_target,
                loss_emphasis=loss_emphases[i],
                n_mirnas=n_mirnas,
                loss_type=loss,
                popt=popt,
                mirnas=miRNAs_filter,
                mirna_expression=mirna_data,
                no_designs=int(designs_per_target/increase_diversity),
                clip_min = clip_min,
                clip_max = clip_max,
                use_combination_penalty=use_combination_penalty,
                forbidden_combinations=forbidden_combinations
            )

            designs["target"] = str(mse_target)
            designs["emphasis"] = str(loss_emphases[i])
            designs["type"] = base_name
        
            used_mirnas = list(designs.index[0])
            single_mirna_knockdown = mirna_data.apply(lambda x: transfer_function_linear(x, *popt))
            knockdown_per_mirna = single_mirna_knockdown.loc[used_mirnas]
            most_potent_mirna = [knockdown_per_mirna.mean(axis=1).idxmin()]

            miRNAs_filter = [mirna for mirna in miRNAs_filter if mirna not in most_potent_mirna]
            target_designs.append(designs.head(int(designs_per_target/increase_diversity)))

        target_designs = pd.concat(target_designs)
        all_designs.append(target_designs)
        
    all_designs_df = pd.concat(all_designs)

    return all_designs_df


def generate_diverse_genetic_designs(base_name, mirna_dataset, targets, emphases, cell_lines_used, popt,
                              loss, sublabel, designs_per_cell_line, n, clip_min=0, clip_max=1,
                            use_combination_penalty=True, forbidden_combinations=None):
    """Generate multiple miRNA circuit designs optimized for different cell type targets,
    iteratively disallowing the most potent miRNA in each generation.
    
    Args:
        base_name (str): Base name used for design naming/indexing
        mirna_dataset (pd.DataFrame): Expression data for miRNAs across cell types
        targets (list): List of target expression patterns to optimize for
        emphases (list): List of dictionaries specifying loss weights for each target
        cell_lines_used (list): Column names of cell types to use from mirna_dataset
        loss (str): Loss function type to use ('mse' or 'binary')
        sublabel (str): Label to add to designs for grouping/filtering
        designs_per_cell_line (int): Number of designs to generate per target
        n (int): Number of miRNAs to use per design
        design_type (str, optional): "active" or "inactive" depending on the desired expression in the cell lines
        use_combination_penalty (bool, optional): Whether to penalize repeated miRNA combinations. Defaults to True.
        forbidden_combinations (list, optional): List of miRNA combinations to avoid. Defaults to None.
    
    Returns:
        pd.DataFrame: Generated circuit designs with performance metrics and metadata
    """
    # we do this here to only get a single design per diversity increase (otherwise, designs tend to be too similar)
    diversity = designs_per_cell_line
    
    designs = generate_multiple_designs(target_values=targets,
                            designs_per_target=designs_per_cell_line, 
                            mirna_data=mirna_dataset,
                            base_name=base_name,
                            loss_emphases=emphases,
                            loss=loss,
                            n_mirnas=n,
                            popt=popt,
                            use_combination_penalty=use_combination_penalty,
                            clip_min=clip_min,
                            clip_max=clip_max,
                            forbidden_combinations=forbidden_combinations,
                            increase_diversity=diversity)
    
    
    designs[cell_lines_used] = designs[cell_lines_used].astype("float")
    designs["sublabel"] = str(sublabel)
    base_number = 4
        
    designs = add_numbered_index(designs, base_name=f"{base_number}_miRNA_{base_name}")
    return designs    


def recalculate_fitnesses(design_df, expression, loss_type="mse", forbidden_combinations=None):
    mirna_columns = [column for column in design_df.columns if "miRNA" in column]
    expression_dataset_with_empty = expression.copy()
    expression_dataset_with_empty.loc["empty"] = 0

    fitness_by_design = {}
    for index, row in design_df.iterrows():
        mirna_tuple = [tuple(row[mirna_columns])]
        emphasis = ast.literal_eval(row["emphasis"])
        mse_target = ast.literal_eval(row["target"])
        target_cell_lines = list(mse_target.keys())
        
        fitness = evaluate_fitness(pop=mirna_tuple, 
                                    expression=expression_dataset_with_empty[target_cell_lines],
                                    loss_emphasis=emphasis,
                                    mse_target=mse_target,
                                    loss_type=loss_type,
                                    forbidden_combinations=forbidden_combinations,
                                    use_combination_penalty=True)
        
        fitness = fitness.reset_index()
        fitness_by_design[index] = fitness
    return pd.concat([fitness_by_design[key].set_index(pd.Index([key])) for key in fitness_by_design.keys()])

def fill_in_other_cell_lines(design_df, expression, popt):
    mirna_columns = [column for column in design_df.columns if "miRNA" in column]
    
    expression_dataset_with_empty = expression.copy()
    expression_dataset_with_empty.loc["empty"] = 0
    
    population = [tuple(row[mirna_columns]) for _, row in design_df.iterrows()]

    # Calculate the stability levels for the designs in the population according to the additive model
    add_expr = add_mirna_combs(expression_dataset_with_empty, population).apply(lambda x: transfer_function_linear(x, *popt))
    design_df[add_expr.columns] = add_expr.values
    
    return design_df