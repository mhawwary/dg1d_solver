


def replace_input_param(fname,str_rep,value_rep):

    n = size(str_rep);

    for lines in fileinput.input(fname, inplace=True):
        a = lines.split('=')
        b = a[0].strip()
        for j in range(0,n):
            if b == str_rep[j]:
                print('{} = {}'.format(a[0].strip(),value_rep[j]))
            else:
                print(lines.rstrip())
