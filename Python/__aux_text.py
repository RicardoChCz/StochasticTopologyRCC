#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:21:40 2019
@author: Ricardo Chávez Cáliz - RCC

Methods to generate .txt information to be used in LaTeX work
"""

def write_table(filename, caption, label, values):
    """
    Write .txt file displaying a table with goodness of fit tests
    Input: filename (str), caption(str), label(str), values(arr)
    Output: .txt file
    """
    # Filename to write
    path = "Txt/"+ filename +".txt"
    # Open the file with writing permission
    myfile = open(path, 'w')
    table_content = ""
    for r in values:
        for c in r:
            table_content +=  '$' + str(c) + '$ &'
        table_content += '\\\\\hline'

    table_begin = """
            \\begin{table}[htbp]
            \\begin{center}
            \\bgroup
            \def\\arraystretch{1.5}
            \\begin{tabular}{|c|c|c|c|}
            \hline
            \diagbox[width=1.3cm, height=0.8cm]{$p$}{$n$} & 5 & 10 & 15 \\\\
            \hline
            """
    table_end = "\end{tabular}\egroup" + "\caption{"+caption+"}"+"\label{"+label+"}"+"\end{center}\end{table}"
    table = table_begin+table_content+table_end
    # Write a line to the file
    myfile.write(table)
 
    # Close the file
    myfile.close()

if __name__ == "__main__":
    array = [
      [1,2,3],
      [3,4,5],
      [6,7,8],
    ]
    write_table("table","caption","label",array)