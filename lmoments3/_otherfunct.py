import scipy as _sp
import scipy.special as _spsp
import scipy.stats as _spst
import math as _math
import sys as _sys

def is_numeric(obj):
    try:
        obj+obj, obj-obj, obj*obj, obj**obj, obj/obj
    except ZeroDivisionError:
        return True
    except Exception:
        return False
    else:
        try:
            a = len(obj)
            if a == 1:
                return True
            return False
        except:
            return True
