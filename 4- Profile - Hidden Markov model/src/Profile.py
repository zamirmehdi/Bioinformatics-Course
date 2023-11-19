import itertools
import math

PSEUDOCOUNT = 2


def calculate_columns(seqs):
    columns = []
    for i in range(len(seqs[0])):
        columns.append([])
        for seq in seqs:
            columns[i].append(seq[i])
    return columns


def calculate_scores():
    for i in range(len(input_seqs[0])):
        # print('\n', i)
        score_matrix[i] = {}
        for char in alphabet:
            char_frequency = seq_columns[i].count(char)
            # print(char_frequency)
            score = (char_frequency + PSEUDOCOUNT) / (len(input_seqs) + len(alphabet) * PSEUDOCOUNT)
            score_matrix[i][char] = score

    char_sum = {}
    for char in alphabet:
        char_sum[char] = 0
    for index in score_matrix:
        for char in score_matrix[index]:
            char_sum[char] += score_matrix[index][char]

    for index in score_matrix:
        for char in alphabet:
            last_score = score_matrix[index][char]
            score_matrix[index][char] = math.log2(last_score / (char_sum.get(char) / 5))


def insert_gaps(word, length, checked_words):
    if len(word) >= length:
        # print([word])
        return [word]
    words = []
    for i in range(len(word) + 1):
        new_word = word[:i] + '-' + word[i:]
        # words.append(new_word)
        if new_word not in checked_words:
            checked_words.append(new_word)
            words += insert_gaps(new_word, length, checked_words)
        # else:
        #     print(word)
    # print(words)
    words = list(dict.fromkeys(words))
    return words


if __name__ == '__main__':
    input_seqs = []
    alphabet = {}
    score_matrix = {}
    number_of_seqs = int(input())

    for i in range(number_of_seqs):
        input_seq = input()
        input_seqs.append(input_seq)

        for char in input_seq:
            if char not in alphabet:
                alphabet[char] = 1
            else:
                alphabet[char] += 1

    search_seq = input()
    seq_columns = calculate_columns(input_seqs)
    calculate_scores()

    # print(input_seqs)
    # print(search_seq)
    # print(alphabet)
    # print(seq_columns)

    cases = []
    max_word = ''
    max_score = -1000
    # for j in range(len(input_seqs[0])-1, -1, -1):
    #     word = search_seq[i:i + j + 1]
    index = 0
    for j in range(len(input_seqs[0]), 2, -1):
        checked_words = []
        for i in range(len(search_seq) - j + 1):
            word = search_seq[i:i + j]
            if word not in checked_words:
                # cases = insert_gaps(word, len(input_seqs[0]), [])
                if len(word) == len(input_seqs[0]):
                    cases = [word]
                else:
                    cases = insert_gaps(word, len(input_seqs[0]), [])
                checked_words.append(word)

                for case in cases:
                    score = 0
                    for k in range(len(case)):
                        score += score_matrix[k][case[k]]
                    if score > max_score:
                        # print(case)
                        max_score = score
                        max_word = case


            # else:
            # print(cases)
            # else:
            #     print(cases)
            #     print(word)

            # for k in range(len(input_seqs[0]) - len(word)):
            #     word += '-'
            # print(word)
            # words = list(map("".join, itertools.combinations(word,len(word))))
            # print(words)
            # cases += words

    # cases = list(dict.fromkeys(cases))
    # print(cases)

    # max_word = ''
    # max_score = -1000
    #
    # for case in cases:
    #     score = 0
    #     for i in range(len(case)):
    #         score += score_matrix[i][case[i]]
    #     if score > max_score:
    #         max_score = score
    #         max_word = case
    print(max_word)

    # print(cases)
    # print()

    # for temp in score_matrix:
    #     print(temp)
    #     print('', score_matrix[temp])
