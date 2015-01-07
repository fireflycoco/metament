import re
__all__ = ['Simple_Expression', 'parse_simple_expression_file']
class Simple_Expression(object):
    """Simple_Expression can deal with expressions with terms in the form of
    'VAR', 'float_coef * VAR' and '+' only.
    e.g.
        A
        B = 0.2 * C + 0.2 * D + E
    """
    def __init__(self, expr_str):
        self.terms = list()
        tokens = self._clean_split(expr_str, sep = '=')
        if len(tokens) == 1:
            self.name = tokens[0]
            self.terms.append((tokens[0], 1.0))
        elif len(tokens) == 2:
            self.name = tokens[0]
            tokens_add = self._clean_split(tokens[1], sep='+')
            for t_add in tokens_add:
                tokens_multi = self._clean_split(t_add, sep='*')
                if len(tokens_multi) == 1:
                    res = self._float_or_str(tokens_multi[0])
                    self.terms.append((res, 1.0))
                elif len(tokens_multi) == 2:
                    res1 = self._float_or_str(tokens_multi[0])
                    res2 = self._float_or_str(tokens_multi[1])
                    self.terms.append((res1, res2))
                elif len(tokens_multi) >= 3:
                    raise ValueError, "Term not acceptable: %s" % tokens_add
        else:
            raise ValueError, "Too many '='s: %s" % expr_str

        match_res = re.search('([^()\s]*)(\s*\(\s*)(.*\S+)(\s*\))', self.name)
        if match_res == None:
            self.unit = ''
        else:
            self.name = match_res.groups()[0]
            self.unit = match_res.groups()[2]
            if len(tokens) == 1:
                self.terms[0] = (self.name, 1.0)

    def __str__(self):
        str_t = [ (str(t1), str(t2)) for t1, t2 in self.terms]
        ts = [' * '.join(t) for t in str_t]
        expr = ' + '.join(ts)
        if self.unit == '':
            return "%s = %s" % (self.name, expr)
        else:
            return "%s ( %s ) = %s" % (self.name, self.unit, expr)

    def eval(self, var_dict):
        result = 0.0
        for t1, t2 in self.terms:
            if type(t1) == str:
                t1 = var_dict[t1]
            if type(t2) == str:
                t2 = var_dict[t2]

            this_term = t1 * t2
            result = result + this_term
        return result

    def var_list(self):
        res = set()
        for term in self.terms:
            for t in term:
                if type(t) == str:
                    res.add(t)
        return list(res)

    def _clean_split(self, s, sep=None):
        tokens = [t.strip() for t in s.split(sep)]
        tokens = [t for t in tokens if t != '']
        return tokens

    def _float_or_str(self, s):
        try:
            res = float(s)
        except ValueError:
            res = s
        return res

def parse_simple_expression_file(fname):
    res = list()
    f = open(fname)
    for line in f:
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue
        expr = Simple_Expression(line)
        res.append(expr)

    return res

