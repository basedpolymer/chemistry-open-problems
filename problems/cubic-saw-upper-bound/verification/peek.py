import struct, sys
f = open(sys.argv[1], 'rb'); data = f.read(40 + 8*200)
print(repr(data[:8]))
d, K = struct.unpack_from('<II', data, 8)
N, P, Q = struct.unpack_from('<QQQ', data, 16)
print('d', d, 'K', K, 'N', N, 'P', P, 'Q', Q)
n = min(N, 30)
codes = struct.unpack_from('<%dQ' % n, data, 40)
for c in codes:
    L = (c >> 57) & 63
    print(hex(c), 'bit63', c >> 63, 'len', L, bin(c & ((1 << 57) - 1)))
