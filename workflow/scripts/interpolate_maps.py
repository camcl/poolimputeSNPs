"""
Script for interpolating (per chromosome) the set of markers in the NIAB Diverse MAGIC wheat population with the IWGSC RefSeq1.0 map.
This script is based on the code in interpolate_maps.py in this repository: https://github.com/joepickrell/1000-genomes-genetic-maps
"""

# TODOs: clean useless parts and document

import os
import numpy as np
import pandas as pd

infilename = snakemake.input.file_id_pos 
mapfilename = snakemake.input.chrom_recomb_map  
tmpoutfilename = snakemake.output.tmp_map 
outfilename = snakemake.output.chrom_interpol_map

posin = list()
varin = list()
mappos = list()
mapgpos = list()

with open(infilename, 'r') as infile: 
        line = infile.readline()
        while line:
                line = line.strip().split('\t')
                pos = int(line[1])
                var = line[0]
                posin.append(pos)
                varin.append(var)
                line = infile.readline()

if os.stat(mapfilename).st_size > 0:
        with open(mapfilename, 'r') as mapfile: 
                line = mapfile.readline()
                while line:
                        line = line.strip().split('\t')
                        pos = int(line[1])
                        gpos = float(line[2])
                        mappos.append(pos)
                        mapgpos.append(gpos)
                        line = mapfile.readline()
else:  # set unit-incremental genetic position (1 cM/Mb = 1e-06/bp as in Beagle4.1) when no markers on the map
        with open(infilename, 'r') as mapfile: 
                line = mapfile.readline()
                gpos = 0.0
                while line:
                        line = line.strip().split('\t')
                        pos = int(line[1])
                        gpos += 1e-06 * (pos - mappos[-1]) if mappos != [] else 0.0
                        mappos.append(pos)
                        mapgpos.append(gpos)
                        line = mapfile.readline()        

with open(tmpoutfilename, 'a') as tmpoutfile:
        index1 = 0
        index2 = 0
        while index1 < len(posin):
                pos = posin[index1]
                var = varin[index1]
                if pos == mappos[index2]:
                        # the MAGIC Wheat site was genotyped as part of the map
                        tmpoutfile.write(f'{var}\t{pos}\t{mapgpos[index2]}\n')
                        index1 = index1+1
                elif pos < mappos[index2]:
                        # current position in interpolation before marker
                        if index2 == 0:
                                # before the first site in the map (genetic position = 0)
                                tmpoutfile.write(f'{var}\t{pos}\t{mapgpos[index2]}\n')
                                index1 = index1+1
                        else:
                                # interpolate
                                prevg = mapgpos[index2-1]
                                prevpos = mappos[index2]
                                frac = (float(pos)-float(mappos[index2-1]))/ (float(mappos[index2]) - float(mappos[index2-1]))
                                tmpg = prevg + frac * (mapgpos[index2] - prevg)
                                tmpoutfile.write(f'{var}\t{pos}\t{tmpg}\n')
                                index1 = index1+1
                elif pos > mappos[index2]:
                        # current position in interpolation after marker
                        if index2 == len(mappos)-1:
                                # after the last site in the map (genetic position = maximum in map, note could try to extrapolate based on rate instead)
                                tmpoutfile.write(f'{var}\t{pos}\t{mapgpos[index2]}\n')
                                index1 = index1+1
                        else:
                                # increment the marker
                                index2 = index2+1

tmpout = pd.read_csv(tmpoutfilename, sep='\t', names=['id', 'pos', 'gpos'])
# gposdiff = tmpout['gpos'].diff().fillna(value=0.0) # 1st gpos is set to 0.0
gposdiff = tmpout['gpos'].apply(lambda pos: pos - tmpout['gpos'].iloc[0])

tmpout['gpos'] = gposdiff
tmpout.to_csv(outfilename, sep=' ', index=False, header=False)
with open(snakemake.log.write_to, 'w') as log_out:
        print(tmpout, file=log_out)

if False:
        ### fix 0s at start of map

        zero_poss = []

        with open(outfilename) as infile:
                line = infile.readline()
                line = line.strip().split()
                dist = float(line[2])
                pos = int(line[1])

                while dist == 0.0:
                        zero_poss.append(pos)
                        line = infile.readline()
                        line = line.strip().split()
                        dist = float(line[2])
                        pos = int(line[1])
                first_nonzero_pos = pos
                first_nonzero_dist = dist


        cM_per_b = first_nonzero_dist/first_nonzero_pos


        outfile = open(outfilename[:-4], 'w')

        with open(outfilename) as infile:
                line = infile.readline()
                outfile.write(line)
                line = line.strip().split()
                dist = float(line[2])
                pos = int(line[1])
                assert dist == 0.0
                for line in infile:
                        line_split = line.strip().split()
                        dist = float(line_split[2])
                        pos = int(line_split[1])
                        if(dist == 0.0):
                                line_split[2] = str(pos*cM_per_b)
                                outfile.write(line_split[0] + " " + line_split[1] + " " + line_split[2] + "\n")
                        else:
                                outfile.write(line)


        os.remove(outfilename)
