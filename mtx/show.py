import fileinput as fi
for line in fi.input():
    for n in line.split():
        print('\033[30m' + n if n == '0' else '\u25fc')
