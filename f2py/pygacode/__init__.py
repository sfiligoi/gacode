import sys

try:
    from gacode_ext import expro

except ImportError:
    __all__ = ['gyro', 'cgyro', 'tgyro', 'neo', 'profiles_gen']

else:
    from gacode_ext import geo
    from gacode_ext import vis

    from . import gyro
    from . import cgyro
    from . import tgyro
    from . import neo
    from . import profiles_gen


    def gapystr_get(s):
        u = []

        try:
            n = len(s)
        except:
            return str(s, 'utf-8')

        for i in range(len(s)):
            if sys.version_info[0] == 2:
                u.append(str(s[i]).strip())
            else:
                u.append(str(s[i], 'utf-8').strip())
        return u


    def gapystr_set(s, l=10, n=20):
        return [item.ljust(l) for item in s + [''] * n][:n]


    __all__ = ['expro', 'geo', 'vis', 'gapystr_get', 'gapystr_set', 'gyro', 'cgyro', 'tgyro', 'neo', 'profiles_gen']
