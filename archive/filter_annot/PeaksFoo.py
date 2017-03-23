import pysam

def annot_intersects(x, v, buf=0):
    annot = []
    for y in v:
        if is_intersect(x, y, buf=buf):
            annot.append(y)
    if len(annot) > 0:
        x.append(annot)
        return x
    else:
        x.append(False)
        return x


def annot_strands(peaks, samfile):
    for k, v in peaks.items():
        new_v = []
        for peak in v:
            start = peak[0]
            end = peak[1]
            strand_counts = get_strands(k, start, end, samfile)
            new_v.append(peak[:-2] + [strand_counts] + peak[-2:])
        peaks[k] = new_v
    return peaks


def filter_annot_peaks(peaks):
    return [peak for peak in peaks if peak[6] or peak[7]]


def filter_peak(p):
    return p[2] > 10 and p[4] > 2


def filter_peaks_dict(d):
    flt_d = {}
    for k, v in d.items():
        flt_d[k] = [p for p in v if filter_peak(p)]
    return flt_d


def get_strands(chrom, start, end, samfile):
    reads = samfile.fetch(chrom, start, end)
    strands = [read.is_reverse for read in reads]
    return  strands.count(False), strands.count(True)


def intersect_peaks(a, b, buf=0):
    c = {}
    for k, v in a.items():
        if k in b:
            c[k] = [peak for peak in v if intersects_known(peak,
                                                           b[k],
                                                           buf=buf)]
    return c

def intersect_peaks_w_extra_annot(a, b, buf=0):
    c = {}
    for k, v in a.items():
        if k in b:
            c[k] = [annot_intersects(peak, b[k], buf=buf) for peak in v]
        else:
            c[k] = [peak + [False] for peak in v]
    return c


def intersects_known(x, v, buf=0):
    for y in v:
        if is_intersect(x, y, buf=buf):
            return True
    return False


def is_intersect(x, y, buf=0):
    return x[1]+buf >= y[0] and x[0]-buf <= y[1]


def merge_dict(a, b):
    c = a.copy()
    for k, v in b.items():
        if k in c:
            c[k] = c[k] + b[k]
        else:
            c[k] = b[k]
    return c


def remove_non_match_peaks(d):
    flt_dict = {}
    for k, v in d.items():
        flt_v = filter_annot_peaks(v)
        if flt_v:
            flt_dict[k] = flt_v
    return flt_dict


def transform_to_dict_w_annot(in_d):
    out_d = {}
    for k, v in in_d.items():
        for peak in v:
            out_d[peak[-4]] = peak
    return out_d


def transform_to_set(d):
    s = set([])
    for k, v in d.items():
        s.update([f[-1] for f in v])
    return s
