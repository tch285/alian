#!/usr/bin/env python
# -*- coding: utf-8 -*-

import alian
import ROOT

mng = ROOT.alian.tracks()
print(mng, type(mng))
ROOT.alian.wait_for_quit()