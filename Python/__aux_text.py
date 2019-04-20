#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 13:21:40 2019
@author: Ricardo Chávez Cáliz - RCC

Methods to generate .txt information to be used in LaTeX work
"""

def write_table():
    """
    Write .txt file displaying a table with goodness of fit tests
    Input: filename (str), caption(str), label(str), values(arr)
    Output: .txt file
    """
    # Filename to write
    filename = "Txt/table.txt"
    # Open the file with writing permission
    myfile = open(filename, 'w')
    
    table = """
            \\begin{table}[htbp]
            \\begin{center}
            \\bgroup
            \def\\arraystretch{1.5}
            \\begin{tabular}{|c|c|c|c|}
            \hline
            \diagbox[width=1.3cm, height=0.8cm]{$p$}{$n$} & 5 & 10 & 15 \\\\
            \hline
             0.2 & $1e^{-8}$ & $1e^{-8}$ & $1e^{-8}$ \\\\\hline
             0.3 & $1e^{-8}$ & $1e^{-8}$ & $1e^{-8}$ \\\\\hline
             0.4 & $1e^{-8}$ & $1e^{-8}$ & $1e^{-8}$ \\\\\hline
            \end{tabular}
            \egroup
            \caption{Maximums of differences between empirical and theoretical probabilities varying $m$ for different values of $n$ and $p$, PYTHON}
            \label{tabla:sencilla}
            \end{center}
            \end{table}
            """

    # Write a line to the file
    myfile.write(table)
 
    # Close the file
    myfile.close()

if __name__ == "__main__":
    write_table()