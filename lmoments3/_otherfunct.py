def is_numeric(obj):
    try:
        obj + obj, obj - obj, obj * obj, obj ** obj, obj / obj
    except ZeroDivisionError:
        return True
    except Exception:
        return False
    else:
        try:
            return bool(len(obj) == 1)
        except:
            return True
