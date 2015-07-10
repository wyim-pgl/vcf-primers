from subprocess import check_output

def boolify(s):
    if s == 'True':
        return True
    if s == 'False':
        return False
    raise ValueError("huh?")

def autoconvert(s):
    for fn in (boolify, int, float):
        try:
            return fn(s)
        except ValueError:
            pass
    return s

def faidx(ref, chrom, start, end):
    loc = "{chrom}:{start}-{end}".format(**locals())
    comm = ["samtools", "faidx", ref, loc]
    return ''.join(check_output(comm).splitlines()[1:])

