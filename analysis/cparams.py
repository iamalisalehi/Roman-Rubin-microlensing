"""Constants of the C++ simulator, read from its own config headers.

The simulator is compiled from include/physical_constants.h, include/common.h, config/data_products.h
and config/parameters.h. This module parses the same text, so the analysis scripts never carry
hand-typed copies of the numbers.

    from cparams import P            # P = load() of this checkout, done at import
    P.Tobs, P.thre[6], P['Rv']       # attribute or item access; arrays come back as numpy arrays
    P.names()                        # every top-level constexpr/const constant found

    python3 analysis/cparams.py            # print every constant, list those it cannot evaluate
    python3 analysis/cparams.py Tobs thre  # print just these

Only top-level `constexpr`/`const` definitions are collected, and each is evaluated lazily (on first
access) by a restricted evaluator: numbers, + - * /, pow/sqrt/float, other constants, lists, and
subscripts. C semantics are kept (int/int truncates; a declared int or double is converted). A
definition it cannot evaluate (struct initialisers such as POPULATIONS) raises CParamError naming the
constant and file:line only when asked for -- there is no silent fallback to a default.
"""
import ast, math, os, re, sys
import numpy as np

HEADERS = ('include/physical_constants.h', 'include/common.h', 'config/data_products.h', 'config/parameters.h')


class CParamError(Exception):
    pass


def _strip(text):
    # string literals -> __S<i>__ placeholders (so '//' or ';' inside them are inert), then comments out
    strs = []
    def keep(m):
        strs.append(ast.literal_eval(m.group(0)))
        return '__S%d__' % (len(strs) - 1)
    text = re.sub(r'"(?:[^"\\\n]|\\.)*"', keep, text)
    text = re.sub(r'/\*.*?\*/', ' ', text, flags=re.S)
    text = re.sub(r'//[^\n]*', '', text)
    text = re.sub(r'^[ \t]*#[^\n]*', '', text, flags=re.M)
    return text, strs


def _split(s):
    # split on commas that are not inside (), {} or []
    out, depth, cur = [], 0, ''
    for ch in s:
        depth += ch in '({[' ; depth -= ch in ')}]'
        if ch == ',' and depth == 0:
            out.append(cur); cur = ''
        else:
            cur += ch
    return out + [cur] if cur.strip() else out


_DECL = re.compile(r'^[ \t]*(?:inline\s+)?(?:static\s+)?(?:constexpr|const)\s+([^;]*);', re.M)
_TYPE = re.compile(r'^((?:const\s+)?(?:std::(?:array|vector)<[^>]*>|[\w:]+)\s*\*?(?:\s*const\b)?)\s*([^=\s].*)$', re.S)
_NUM = re.compile(r'(\b\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?)[fFlLuU]+\b')


class Params:
    def __init__(self, root=None):
        self.root = root or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self._defs, self._val, self._busy = {}, {}, set()
        for h in HEADERS:
            self._read(h)

    def _read(self, rel):
        raw = open(os.path.join(self.root, rel)).read()
        text, self._strs = _strip(raw)
        depth, prev = 0, 0
        for m in _DECL.finditer(text):
            seg = text[prev:m.start()]
            depth += seg.count('{') - seg.count('}'); prev = m.start()
            if depth != 0:
                continue                      # a local inside a function body
            line = text.count('\n', 0, m.start()) + 1
            t = _TYPE.match(m.group(1).strip())
            if not t:
                continue
            ctype = re.sub(r'\s+', '', re.sub(r'^const\s+', '', t.group(1)))   # 'int', 'char*', 'char*const', 'std::array<double,M>'
            for d in _split(t.group(2)):
                n = re.match(r'\s*(\w+)\s*(\[[^\]]*\])?\s*=\s*(.*)$', d, re.S)
                if n:
                    expr = re.sub(r'__S(\d+)__', lambda k: repr(self._strs[int(k.group(1))]), n.group(3).strip())
                    self._defs[n.group(1)] = (ctype, bool(n.group(2)), expr, '%s:%d' % (rel, line))

    def _eval(self, name, node):
        ev = lambda n: self._eval(name, n)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float, str)):
            return node.value
        if isinstance(node, ast.List):
            return [ev(e) for e in node.elts]
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            v = ev(node.operand)
            return -v if isinstance(node.op, ast.USub) else v
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div)):
            a, b = ev(node.left), ev(node.right)
            if isinstance(node.op, ast.Add): return a + b
            if isinstance(node.op, ast.Sub): return a - b
            if isinstance(node.op, ast.Mult): return a * b
            if isinstance(a, int) and isinstance(b, int):
                return (abs(a) // abs(b)) * (1 if (a < 0) == (b < 0) else -1)   # truncate toward zero
            return a / b
        if isinstance(node, ast.Name):
            return math.pi if node.id == 'M_PI' else self[node.id]
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and not node.keywords:
            f = {'pow': math.pow, 'float': float, 'sqrt': math.sqrt}.get(node.func.id)
            if f:
                return f(*[ev(a) for a in node.args])
        if isinstance(node, ast.Subscript):
            return ev(node.value)[int(ev(node.slice))]
        raise CParamError('%s (%s): cannot evaluate %s' % (name, self._defs[name][3], ast.unparse(node)))

    def _value(self, name):
        ctype, is_arr, expr, where = self._defs[name]
        e = re.sub(r'\bstd::', '', expr)
        e = re.sub(r'\bdouble\s*\(', 'float(', e)
        e = _NUM.sub(r'\1', e).replace('{', '[').replace('}', ']')
        try:
            v = self._eval(name, ast.parse(e, mode='eval').body)
        except SyntaxError:
            raise CParamError('%s (%s): not an arithmetic expression: %s' % (name, where, expr[:60]))
        if is_arr or ctype.startswith('std::'):
            v = v if isinstance(v, list) else [v]
            if ctype.startswith('std::array<double'):
                return np.array([float(x) for x in v])
            if ctype.startswith('std::vector<int'):
                return [int(x) for x in v]
            return tuple(v)
        return {'int': int, 'double': float}.get(ctype, lambda x: x)(v)

    def __getitem__(self, name):
        if name not in self._defs:
            raise AttributeError('%s is not defined in config/*.h, include/common.h or include/physical_constants.h' % name)
        if name not in self._val:
            if name in self._busy:
                raise CParamError('%s (%s): circular definition' % (name, self._defs[name][3]))
            self._busy.add(name)
            try:
                self._val[name] = self._value(name)
            finally:
                self._busy.discard(name)
        return self._val[name]

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        return self[name]

    def names(self):
        return list(self._defs)


def load(root=None):
    return Params(root)


P = load()

if __name__ == '__main__':
    todo = sys.argv[1:] or P.names()
    bad = []
    for n in todo:
        try:
            print('%-28s %s' % (n, P[n] if not isinstance(P[n], np.ndarray) else P[n].tolist()))
        except (CParamError, AttributeError) as e:
            bad.append(str(e))
    if bad:
        print('\ncannot evaluate:', *bad, sep='\n  ')
