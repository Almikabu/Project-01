import string
from time import time

import utils
from utils import *


rounds_no = None

chromosome_length = None

decryption_key = None


def set_globals():
    global rounds_no
    global decryption_key
    rounds_no = random.randrange(3, 12, 2)
    decryption_key = ""


def encrypt_key(data, key):
    if len(data) > len(key):
        factor = int(len(data) / len(key))
        key += key * factor

        return bitxor(data, key)

    return bitxor(data, key)


def reshape(dna_sequence):
    
    global chromosome_length
    global decryption_key

    divs = divisors(len(dna_sequence))
    chromosome_no = divs[random.randint(0, len(divs) - 1)]
    chromosome_length = int(len(dna_sequence) / chromosome_no)
    chromosomes = []

    decryption_key += reshape_del + str(chromosome_length) + reshape_del

    for i in range(0, len(dna_sequence), chromosome_length):
        chromosomes.append(dna_sequence[i:i + chromosome_length])

    return chromosomes


def reverse_reshape(population):
    return "".join(population)


def rotate_crossover(population):
    
    global chromosome_length
    global decryption_key

    new_population = []

    decryption_key += rotate_crossover_del

    rotation_offset = random.randint(1, chromosome_length)

    decryption_key += rotation_offset_del + str(rotation_offset) + rotation_offset_del

    decryption_key += rotation_types_del

    for chromosome in population:

        p = random.uniform(0, 1)

        if p > 0.5:
            decryption_key += "right|"
            right_first = chromosome[0: len(chromosome) - rotation_offset]
            right_second = chromosome[len(chromosome) - rotation_offset:]
            new_population.append(right_second + right_first)
        else:
            decryption_key += "left|"
            left_first = chromosome[0: rotation_offset]
            left_second = chromosome[rotation_offset:]
            new_population.append(left_second + left_first)

    decryption_key += rotation_types_del

    decryption_key += rotate_crossover_del

    return new_population


def single_point_crossover(population):
    
    global decryption_key

    decryption_key += single_point_crossover_del

    new_population = []
    for i in range(0, len(population) - 1, 2):
        candidate1 = population[i]
        candidate2 = population[i + 1]

        length = len(candidate1)
        crossover_point = random.randint(0, length - 1)

        decryption_key += str(crossover_point) + "|"

        offspring1 = candidate2[0: crossover_point] + candidate1[crossover_point:]
        offspring2 = candidate1[0: crossover_point] + candidate2[crossover_point:]
        new_population.append(offspring1)
        new_population.append(offspring2)

    if len(population) % 2 == 1:
        new_population.append(population[len(population) - 1])

    decryption_key += single_point_crossover_del

    return new_population


def crossover(population):
    global decryption_key

    p = random.uniform(0, 1)

    if p < 0.33:
        decryption_key += crossover_type_del + "rotate_crossover" + crossover_type_del
        return rotate_crossover (population)
    elif p >= 0.33 and p < 0.66:
        decryption_key += crossover_type_del + "single_point_crossover" + crossover_type_del
        return single_point_crossover(population)
    else:
        decryption_key += crossover_type_del + "both" + crossover_type_del
        population = rotate_crossover (population)
        return single_point_crossover(population)


def complement(chromosome, point1, point2):
    
    new_chromosome = ""

    for i in range(len(chromosome)):
        if i >= point1 and i <= point2:
            if chromosome[i] == '0':
                new_chromosome += '1'
            else:
                new_chromosome += '0'
        else:
            new_chromosome += chromosome[i]

    return new_chromosome


def alter_dna_bases(bases):
    
    alter_dna_table = {}

    for _ in range(2):
        base1 = bases[random.randint(0, len(bases) - 1)]
        bases.remove(base1)

        base2 = bases[random.randint(0, len(bases) - 1)]
        bases.remove(base2)

        alter_dna_table[base1] = base2
        alter_dna_table[base2] = base1

    return alter_dna_table


def mutation(population):
    
    global decryption_key

    bases = ['A', 'C', 'G', 'T']
    alter_dna_table = alter_dna_bases(bases)

    decryption_key += mutation_table_del + str(alter_dna_table) + mutation_table_del

    new_population = []
    for chromosome in population:
        decryption_key += chromosome_del

        b_chromosome = dna_to_bits(chromosome, utils.dna_base_to_two_bits_table)
        decryption_key += complement_mutation_del
        point1 = random.randint(0, len(b_chromosome) - 1)
        point2 = random.randint(point1, len(b_chromosome) - 1)
        decryption_key += "(%s, %s)" % (point1, point2)
        decryption_key += complement_mutation_del
        b_chromosome = complement(b_chromosome, point1, point2)

        four_bits_vector = group_bits(b_chromosome, 4)

        last_dna_base = None
        
        if len(four_bits_vector[len(four_bits_vector) - 1]) == 2:
            last_dna_base = utils.two_bits_to_dna_base_table[four_bits_vector[len(four_bits_vector) - 1]]

            four_bits_vector = four_bits_vector[:-1]

        dna_seq = bits_to_dna(four_bits_vector, utils.four_bits_to_two_dna_base_table)
        if last_dna_base is not None:
            dna_seq += last_dna_base

        decryption_key += alter_mutation_del
        point1 = random.randint(0, len(dna_seq) - 1)
        point2 = random.randint(point1, len(dna_seq) - 1)
        decryption_key += "(%s, %s)" % (point1, point2)
        decryption_key += alter_mutation_del
        new_chromosome = ""
        for i in range(len(dna_seq)):
            if i >= point1 and i <= point2:
                new_chromosome += alter_dna_table[dna_seq[i]]
            else:
                new_chromosome += dna_seq[i]

        new_population.append(new_chromosome)

        decryption_key += chromosome_del

    return new_population


def dnaEncrypt(text, key):
    global rounds_no
    global decryption_key
    global decryption_key
    
    
    decryption_key = key_del + key + key_del

    print("\nDNA Encryption is running...\n")

    b_data1 = binarized_data(text)
    dna_seq = bits_to_dna(b_data1, utils.two_bits_to_dna_base_table)

    b_data2 = dna_seq
    print("Initial DNA sequence:", dna_seq)

    decryption_key += no_rounds_del + str(rounds_no) + no_rounds_del

    while rounds_no > 0:
        decryption_key += round_del

        b_data2 = bits_to_dna(
            group_bits(encrypt_key(dna_to_bits(reverse_reshape(b_data2), utils.dna_base_to_two_bits_table), key)),
            utils.two_bits_to_dna_base_table)
        
        b_data2 = reshape(b_data2)
        

        decryption_key += crossover_del
        b_data2 = crossover(b_data2)
        decryption_key += crossover_del
        
        decryption_key += mutation_del
        b_data2 = mutation(b_data2)
        decryption_key += mutation_del

        rounds_no -= 1

        decryption_key += round_del

    return reverse_reshape(b_data2)


def main():
    global decryption_key

    original_file = open(original_filename, "w")

    text = "I want to eat pinnaple with satria indah"

    print("Text:", text)

    original_file.write(text)

    key = str2bin(''.join(random.SystemRandom().choice(string.ascii_letters + string.digits) for _ in range(16)))

    print("Key:", len(key), key)

    set_globals()

    decryption_key += key_del + key + key_del

    generate_pre_processing_tables()
    generate_mutation_tables()

    start = time()
    encrypted_text = dnaEncrypt(text, key)
    print("Final DNA sequence:", encrypted_text)
    end = time()

    print("\nTotal execution time:", end - start)

    key_file = open(key_filename, "w")
    encrypted_file = open(encrypted_filename, "w")

    encrypted_file.write(encrypted_text)

    key_file.write(decryption_key)

    encrypted_file.close()
    original_file.close()
    key_file.close()


if __name__ == '__main__':
    main()
