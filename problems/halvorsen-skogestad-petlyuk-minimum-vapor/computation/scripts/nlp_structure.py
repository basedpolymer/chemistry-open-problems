"""Configuration structure helper (no solver dependency).

structure(cfg) returns (cols, info, reb, side):
  cols: mixture streams (columns); info[s]: dict(i, j, k, l, top, bot, kind) where [i,k]/[l,j] is the split,
  top/bot the producing column (distillate/residue link) or None, kind in {'TC','HX','SD'};
  reb: list of columns whose bottom vapor is generated in a reboiler (submixture reboilers and pure products
  produced only as residues); side: list of (top producer, bottom producer) of pure side-draw products.
Environment: Python 3 standard library.
Runtime: negligible.
"""


def structure(cfg):
    """Precompute per-column data."""
    N = cfg['N']
    cols = list(cfg['streams'])
    info = {}
    for s in cols:
        i, j = s
        k, l = cfg['splits'][s]
        prods = cfg['producers'].get(s, [])
        top = [p for (t, p) in prods if t == 'top']
        bot = [p for (t, p) in prods if t == 'bot']
        kind = cfg['hx'].get(s, 'SD')  # 'TC', 'HX' or 'SD' (side draw) ; [1,N] -> 'SD' unused
        info[s] = dict(i=i, j=j, k=k, l=l, top=top[0] if top else None, bot=bot[0] if bot else None, kind=kind)
    # reboilers: (bottom producer column) list
    reb = []
    for s in cols:
        if s == (1, N):
            continue
        d = info[s]
        if d['top'] is None and d['kind'] == 'HX':
            reb.append(d['bot'])
    side = []
    for p in range(1, N + 1):
        prods = cfg['producers'][(p, p)]
        top = [x for (t, x) in prods if t == 'top']
        bot = [x for (t, x) in prods if t == 'bot']
        if top and bot:
            side.append((top[0], bot[0]))
        elif bot:
            reb.append(bot[0])
    return cols, info, reb, side
