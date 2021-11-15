import os
math_ops = ["*", "+", "-", "/"]


def next_val(ar, x):
    for i, v in enumerate(ar):
        if v > x:
            return i
    raise ValueError


def prev_val(ar, x):
    for i, v in enumerate(ar[::-1]):
        if v < x:
            return i
    raise ValueError

def mat_to_c(line):
    pow_indices = [i for i, x in enumerate(line) if x == "^"]
    math_indices = [i for i, x in enumerate(line) if x in math_ops]

    ob_indices = [i for i, x in enumerate(line) if x == "("]
    cb_indices = [i for i, x in enumerate(line) if x == ")"]
    assert len(ob_indices) == len(cb_indices)

    for pi in pow_indices:
        # ()^
        if pi - 1 in cb_indices:
            start = prev_val(ob_indices, pi-1)
            lhs = line[start + 1: pi-1]

            # ()^()
            if pi + 1 in ob_indices:
                end = next_val(cb_indices, pi+1)
                rhs = line[pi+2: end]
            # ()^M
            else:
                end = next_val(math_indices, pi) - 1
                rhs = line[pi+1: end+1]
        # M^
        else:
            start = prev_val(math_indices, pi) + 1
            lhs = line[start:pi]

            # M^()
            if pi + 1 in ob_indices:
                end = next_val(cb_indices, pi+1)
                rhs = line[pi+2: end]

            # M^M
            else:
                end = next_val(math_indices, pi) - 1
                rhs = line[pi+1: end+1]

        print(f"String processed: {line[start:end+1]}")
        print(f"LHS: {lhs} \t RHS: {rhs}")
        import pdb
        pdb.set_trace()


file_path = "d2psi_dF2.txt"
lines = open(file_path).readlines()
# lines = [l.replace(" ", "") for l in lines] # remove spaces
lines = [l.strip() for l in lines] # remove spaces
c_lines = [mat_to_c(l) for l in lines]

