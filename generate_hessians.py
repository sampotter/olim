import sympy

tetras = [
    (int('011', 2), int('010', 2), int('110', 2)),
    (int('110', 2), int('100', 2), int('101', 2)),
    (int('101', 2), int('001', 2), int('011', 2)),
    (int('001', 2), int('011', 2), int('110', 2)),
    (int('011', 2), int('010', 2), int('100', 2)),
    (int('010', 2), int('110', 2), int('101', 2)),
    (int('110', 2), int('100', 2), int('001', 2)),
    (int('100', 2), int('101', 2), int('011', 2)),
    (int('101', 2), int('001', 2), int('010', 2)),
    (int('001', 2), int('011', 2), int('100', 2)),
    (int('011', 2), int('010', 2), int('101', 2)),
    (int('010', 2), int('110', 2), int('001', 2)),
    (int('110', 2), int('100', 2), int('011', 2)),
    (int('100', 2), int('101', 2), int('010', 2)),
    (int('101', 2), int('001', 2), int('110', 2)),
    (int('001', 2), int('010', 2), int('100', 2)),
    (int('011', 2), int('110', 2), int('101', 2)),
    (int('001', 2), int('011', 2), int('111', 2)),
    (int('011', 2), int('010', 2), int('111', 2)),
    (int('010', 2), int('110', 2), int('111', 2)),
    (int('110', 2), int('100', 2), int('111', 2)),
    (int('100', 2), int('101', 2), int('111', 2)),
    (int('101', 2), int('001', 2), int('111', 2)),
    (int('001', 2), int('010', 2), int('111', 2)),
    (int('010', 2), int('100', 2), int('111', 2)),
    (int('100', 2), int('001', 2), int('111', 2)),
    (int('011', 2), int('110', 2), int('111', 2)),
    (int('110', 2), int('101', 2), int('111', 2)),
    (int('101', 2), int('011', 2), int('111', 2)),
    # these are for update_rules.tetra_updates.test.cpp
    (int('001', 2), int('011', 2), int('101', 2)),
    (int('001', 2), int('101', 2), int('110', 2)),
    (int('011', 2), int('101', 2), int('110', 2))]

def rational_cstr(r):
    if r == 0:
        return '0'
    elif r.is_integer:
        return str(r)
    else:
        return '%.1f/%.1f' % r.as_numer_denom()

def get_bit(i, pos):
    return (i & (1 << pos)) >> pos

fmt = 'template <> int %s<%d, %d, %d, 2>::_dPt_dP[] = {%s, %s, %s};'

cstrs = dict()

for i, j, k in tetras:
    print((i, j, k))
    P = sympy.Matrix([
        [get_bit(i, 0), get_bit(j, 0), get_bit(k, 0)],
        [get_bit(i, 1), get_bit(j, 1), get_bit(k, 1)],
        [get_bit(i, 2), get_bit(j, 2), get_bit(k, 2)]])
    dP = P[:, 1:] - P[:, 0].row_join(P[:, 0])
    dPt_dP = dP.T*dP
    a_cstr = rational_cstr(dPt_dP[0, 0])
    b_cstr = rational_cstr(dPt_dP[0, 1])
    c_cstr = rational_cstr(dPt_dP[1, 1])
    cstrs[(i, j, k)] = (a_cstr, b_cstr, c_cstr)

for i, j, k in tetras:
    a_cstr, b_cstr, c_cstr = cstrs[(i, j, k)]
    print(fmt % ('F0_bv', i, j, k, a_cstr, b_cstr, c_cstr))

print()

for i, j, k in tetras:
    a_cstr, b_cstr, c_cstr = cstrs[(i, j, k)]
    print(fmt % ('F1_bv', i, j, k, a_cstr, b_cstr, c_cstr))
