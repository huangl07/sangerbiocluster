# -*- coding:utf-8 -*-
"""
LastEditTime: 2023/06/06
Author: yuan.xu
mail: yuan.xu@majorbio.com
"""

import os
import argparse
from subprocess import Popen


def run_cmd(shell_script,text):
    shell_script = f"sh {shell_script}"
    subp = Popen(shell_script,shell=True,encoding="utf-8")
    subp.wait()
    if subp.poll() == 0:
        print("{}: Run Success!".format(text))
    else:
        print("{}: Run Error!".format(text))
