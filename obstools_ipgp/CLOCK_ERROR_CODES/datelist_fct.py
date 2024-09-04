from obspy.core import UTCDateTime


def datelist(start, end, delta):
    """
    usage: 
    from datetime import timedelta
    datevec_new = perdelta(datevec[0], datevec[-1], timedelta(days=1)
    """
    datelist = []
    curr = start
    while curr <= end:
        curr2 = round(curr, 3)
        curr = UTCDateTime(curr2)
        datelist.append(curr)
        curr += delta
    return datelist