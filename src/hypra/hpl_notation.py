"""Exact HPL index and plus/minus kernel notation conversions."""
from itertools import product

import sympy as s

from .errors import UnsupportedExpressionError
from .hpl_algebra import _word
from .symbolic import HPL


def expand_indices(indices):
    """Expand signed compressed indices; a standalone zero stays a zero."""
    if not isinstance(indices,(tuple,list,s.Tuple)) or any(
        isinstance(n,bool) or not isinstance(n,(int,s.Integer)) for n in indices
    ):
        raise TypeError("Compressed indices must be exact integers")
    if sum(max(1,abs(int(n))) for n in indices)>24:
        raise UnsupportedExpressionError("Index conversion weight exceeds24")
    return tuple(letter for n in indices
                 for letter in ((0,) if n==0 else (0,)*(abs(int(n))-1)+(int(s.sign(n)),)))


def compress_indices(word):
    """Canonical compression: zero runs attach to the next nonzero letter."""
    word=_word(word)
    result,zeros=[],0
    for letter in word:
        if letter==0:
            zeros+=1
        else:
            result.append(letter*(zeros+1))
            zeros=0
    return tuple(result)+(0,)*zeros


def _linear_choices(choices):
    if len(choices)>12:
        raise UnsupportedExpressionError("Plus/minus conversion weight exceeds12")
    result={}
    for terms in product(*choices):
        word=tuple(letter for letter,_ in terms)
        result[word]=s.prod(coefficient for _,coefficient in terms)
    return result


def pm_to_letters(word):
    """Word coefficient dictionary for f_plus=f_1+f_-1, f_minus=f_1-f_-1."""
    if not isinstance(word,(tuple,list,s.Tuple)):
        raise TypeError("Plus/minus word must be a sequence")
    choices=[]
    for letter in word:
        if type(letter) is int and letter==0 or isinstance(letter,s.Integer) and letter==0:
            choices.append(((0,s.S.One),))
        elif letter=='plus':
            choices.append(((1,s.S.One),(-1,s.S.One)))
        elif letter=='minus':
            choices.append(((1,s.S.One),(-1,-s.S.One)))
        else:
            raise ValueError("Use zero, 'plus', or 'minus' as kernel letters")
    return _linear_choices(choices)


def letters_to_pm(word):
    """Inverse word dictionary: f_±1=(f_plus ± f_minus)/2."""
    choices=[]
    for letter in _word(word):
        choices.append(((0,s.S.One),) if letter==0 else
                       (('plus',s.S.Half),('minus',s.Rational(letter,2))))
    return _linear_choices(choices)


def hpl_pm(word, argument):
    """Construct plus/minus notation as an exact sum of ordinary HPL objects."""
    if not isinstance(argument,s.Expr):
        raise TypeError("Argument must be a SymPy expression")
    return sum((coefficient*(HPL(letters,argument) if letters else s.S.One)
                for letters,coefficient in pm_to_letters(word).items()),s.S.Zero)
