import os


def save_as_multiple_file(number, size, field, file_type, fieldtype):
    idx = number.strip('0')
    os.makedirs(f'./{fieldtype}', exist_ok=True)
    idx = 0 if not idx else int(idx)
    with open(f'./{fieldtype}/{fieldtype}{number}.{file_type}', 'w') as f:
        for j in range(0, size):
            print(
                '%.2d' % (j),
                '%10.6f' % (field[idx, j]),
                file=f,
            )


def save_as_one_file(path, field):
    with open(f'{path}', 'w') as f:
        for i in range(len(field[:, 0])):
            for j in range(len(field[0, :])):
                end = '\n' if j == len(field[0, :]) - 1 else ', '
                print(
                    '%10.6f' % (field[i, j]),
                    file=f,
                    end=end,
                )
