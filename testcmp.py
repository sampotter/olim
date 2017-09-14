#!/usr/bin/env python3

import fileinput

def is_failure_line(line):
    return line[:7] == 'failure'

def is_segfault_line(line):
    return line[:18] == 'segmentation fault'

def process_line(line):
    if is_failure_line(line) or is_sefault_line(line):
        sys.exit("test failure")

if __name__ == '__main__':
    map(process_line(line) for line in fileinput.input())
        
