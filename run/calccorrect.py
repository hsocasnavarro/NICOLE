#!/usr/bin/env python
# -*- coding: utf-8 -*-

def message(correct=True):
    import curses

    stdscr = curses.initscr()
    [maxy,maxx]=stdscr.getmaxyx()
    curses.endwin()

    line1=' CALCULATIONS '
    if correct:
        line2='   CORRECT    '
    else:
        line2='  INCORRECT!  '
    len1=len(line1)
    len2=len(line2)
    maxlen=max(len1,len2)
    sizey=4
    sizex=maxlen+2
    y0=int(maxy/2.)-2
    x0=int(maxx/2.)-int(sizex/2.)
    sx0=str(x0)
    sy0=str(y0)

    shorizontalup=u'\u250F'+u'\u2501'*(maxlen)+u'\u2513'
    sline1=u'\u2503'+line1+u'\u2503'
    sline2=u'\u2503'+line2+u'\u2503'
    shorizontaldown=u'\u2517'+u'\u2501'*(maxlen)+u'\u251B'


    print("\033["+sy0+";"+sx0+"H"+shorizontalup).encode('utf-8')
    sy0=str(y0+1)
    print("\033["+sy0+";"+sx0+"H"+sline1).encode('utf-8')
    sy0=str(y0+2)
    print("\033["+sy0+";"+sx0+"H"+sline2).encode('utf-8')
    sy0=str(y0+3)
    print("\033["+sy0+";"+sx0+"H"+shorizontaldown).encode('utf-8')

    # Position cursor
    sy0=str(maxy-1)
    sx0=str(1)
    print("\033["+sy0+";"+sx0+"H").encode('utf-8')

