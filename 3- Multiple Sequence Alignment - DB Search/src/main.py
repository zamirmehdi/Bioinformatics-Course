import operator

S_MATCH = 3
S_MISSMATCH = -1
S_GAP = -2


def global_align(x, y, s_match, s_mismatch, s_gap):
    A = []

    for i in range(len(y) + 1):
        A.append([0] * (len(x) + 1))

    for i in range(len(y) + 1):
        A[i][0] = s_gap * i

    for i in range(len(x) + 1):
        A[0][i] = s_gap * i

    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            A[i][j] = max(
                A[i][j - 1] + s_gap,
                A[i - 1][j] + s_gap,
                A[i - 1][j - 1] + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != '-') else 0) + (
                    s_mismatch if (y[i - 1] != x[j - 1] and y[i - 1] != '-' and x[j - 1] != '-') else 0) + (
                    s_gap if (y[i - 1] == '-' or x[j - 1] == '-') else 0)
            )
    align_X = ""
    align_Y = ""

    i = len(x)
    j = len(y)

    while i > 0 or j > 0:
        current_score = A[j][i]

        if i > 0 and j > 0 and (
                ((x[i - 1] == y[j - 1] and y[j - 1] != '-') and current_score == A[j - 1][i - 1] + s_match) or
                ((y[j - 1] != x[i - 1] and y[j - 1] != '-' and x[i - 1] != '-') and current_score == A[j - 1][
                    i - 1] + s_mismatch) or
                ((y[j - 1] == '-' or x[i - 1] == '-') and current_score == A[j - 1][i - 1] + s_gap)
        ):

            align_X = x[i - 1] + align_X
            align_Y = y[j - 1] + align_Y
            i = i - 1
            j = j - 1

        elif i > 0 and (current_score == A[j][i - 1] + s_gap):
            align_X = x[i - 1] + align_X
            align_Y = "-" + align_Y
            i = i - 1

        else:
            align_X = "-" + align_X
            align_Y = y[j - 1] + align_Y
            j = j - 1

    return (align_X, align_Y, A[len(y)][len(x)])


def get_input():
    number_of_seqs = input()
    input_seqs = []
    for i in range(int(number_of_seqs)):
        input_seqs.append(input())
    return input_seqs


def fill_matrix_and_find_center(seqs):
    max_score = float('-inf')
    center_seq = ''
    score_matrix = {}

    for seq1 in seqs:
        score_matrix[seqs.index(seq1)] = {}
        total_score = 0

        for seq2 in seqs:

            if seqs.index(seq1) != seqs.index(seq2):
                x, y, score = global_align(seq1, seq2, S_MATCH, S_MISSMATCH, S_GAP)
                score_matrix[seqs.index(seq1)][seqs.index(seq2)] = score

                total_score += score

        if total_score > max_score:
            max_score = total_score

            center_seq = seq1

    return score_matrix, center_seq


def always_a_gap(center, updated_center, last_seqs):
    loc1 = 0
    loc2 = 0
    loc3 = 0
    total_center = ''

    while loc1 < len(center) and loc2 < len(updated_center):

        if center[loc1] == updated_center[loc2]:
            total_center += center[loc1]
            loc1 += 1
            loc2 += 1
            loc3 += 1

        elif center[loc1] == '-':
            total_center += '-'
            loc1 += 1
            loc3 += 1

        elif updated_center[loc2] == '-':
            total_center += '-'
            for seq in last_seqs:
                for index in range(len(seq)):
                    if index == loc3:
                        last_seqs[last_seqs.index(seq)] = seq[0:loc3] + '-' + seq[loc3:]
                        break
            loc2 += 1
            loc3 += 1

    for i in range(loc3, loc3 + len(updated_center) - loc2):
        for seq in last_seqs:
            total_center += '-'
            last_seqs[last_seqs.index(seq)] = seq[0:i] + '-' + seq[i:]

    for i in range(loc3, loc3 + len(updated_center) - loc1):
        for seq in last_seqs:
            total_center += '-'


def star_alignment(center_seq, score_matrix, seqs):
    sorted_seqs = sorted(score_matrix[seqs.index(center_seq)].items(), key=operator.itemgetter(1), reverse=True)

    aligned_seqs = []
    new_center_seq = center_seq

    for seq in sorted_seqs:
        last_center_seq = new_center_seq

        new_seq, new_center_seq, score = global_align(seqs[seq[0]], new_center_seq, S_MATCH, S_MISSMATCH, S_GAP)

        always_a_gap(last_center_seq, new_center_seq, aligned_seqs)
        aligned_seqs.append(new_seq)

    aligned_seqs = [new_center_seq] + aligned_seqs

    output_seqs = []
    for seq in seqs:
        for aligned_seq in aligned_seqs:
            if aligned_seq.replace('-', '') == seq:
                output_seqs.append(aligned_seq)

    return output_seqs


def calculate_alignment_score(columns):
    total_score = 0
    for column in columns:
        for i in range(len(column)):
            for j in range(i + 1, len(column)):
                if column[i] == '-' and column[j] == '-':
                    total_score += 0
                elif column[i] == '-' or column[j] == '-':
                    total_score += -2
                elif column[i] != column[j]:
                    total_score += -1
                elif column[i] == column[j]:
                    total_score += 3
    return total_score


def calculate_columns(aligned_seqs):
    columns = []
    for i in range(len(aligned_seqs[0])):
        columns.append([])
        for aligned_seq in aligned_seqs:
            columns[i].append(aligned_seq[i])
    return columns


def block_improvement(seqs_columns, aligned_seqs):
    improved = False
    not_valid_column_indexes = []
    for i in range(len(seqs_columns)):
        same_chars = 0
        for j in range(len(seqs_columns[i])):
            if seqs_columns[i][j] != seqs_columns[i][j - 1]:
                same_chars = 1
                break
        if same_chars == 0:
            not_valid_column_indexes.append(i)

    blocks = []
    for i in range(len(not_valid_column_indexes)):

        if i == len(not_valid_column_indexes) - 1:
            if not_valid_column_indexes[i] - not_valid_column_indexes[i - 1] > 2:
                blocks.append([not_valid_column_indexes[i - 1] + 1, not_valid_column_indexes[i]])
                # print([not_valid_column_indexes[i - 1] + 1, not_valid_column_indexes[i]])
            if len(seqs_columns) - 1 - not_valid_column_indexes[i] > 2:
                blocks.append([not_valid_column_indexes[i] + 1, len(seqs_columns)])
            # print([not_valid_column_indexes[i] + 1, len(seqs_columns)])

        elif i == 0:
            if not_valid_column_indexes[i] > 1:
                blocks.append([0, not_valid_column_indexes[i]])
                # print([0, not_valid_column_indexes[i] - 1])

        elif not_valid_column_indexes[i] - not_valid_column_indexes[i - 1] > 2:
            blocks.append([not_valid_column_indexes[i - 1] + 1, not_valid_column_indexes[i]])
            # print([not_valid_column_indexes[i - 1] + 1, not_valid_column_indexes[i]])

    final_seqs = aligned_seqs.copy()
    len_difference = 0

    for block in blocks:
        # print(block)
        # print(seqs_columns[block[0]:block[1]])
        last_score = calculate_alignment_score(columns=seqs_columns[block[0]:block[1]])
        # print(last_score)

        block_seqs = []
        for seq in aligned_seqs:
            block_seqs.append(seq[block[0]:block[1]].replace('-', ''))

        block_score_matrix, block_center_seq = fill_matrix_and_find_center(block_seqs)
        aligned_block_seqs = star_alignment(block_center_seq, block_score_matrix, block_seqs)
        block_seqs_columns = calculate_columns(aligned_block_seqs)
        block_alignment_score = calculate_alignment_score(block_seqs_columns)
        # print(block_alignment_score)

        if block_alignment_score > last_score:
            improved = True
            new_alignment = []
            for i in range(len(final_seqs)):
                final_seqs[i] = final_seqs[i][0:block[0] - len_difference] + aligned_block_seqs[i] + final_seqs[i][
                                                                                                     block[
                                                                                                         1] - len_difference:]
                # print(final_seqs)

            len_difference = block[1] - block[0] - len(aligned_block_seqs[0])
    return improved, final_seqs


if __name__ == '__main__':
    seqs = get_input()
    updated_center_seq = ''
    score_matrix, center_seq = fill_matrix_and_find_center(seqs)

    aligned_seqs = star_alignment(center_seq, score_matrix, seqs)
    seqs_columns = calculate_columns(aligned_seqs)
    alignment_score = calculate_alignment_score(seqs_columns)

    print(alignment_score)
    for output_seq in aligned_seqs:
        print(output_seq)
    print()

    output_seqs = aligned_seqs

    while True:
        # print()
        output_columns = calculate_columns(output_seqs)
        # print('aaa\n', output_seqs)
        improved, output_seqs = block_improvement(output_columns, output_seqs)
        if not improved:
            break

    output_columns = calculate_columns(output_seqs)
    improved, output_seqs = block_improvement(output_columns, output_seqs)

    output_columns = calculate_columns(output_seqs)
    output_score = calculate_alignment_score(output_columns)
    print(output_score)
    for output_seq in output_seqs:
        print(output_seq)
