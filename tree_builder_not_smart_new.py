from ete3 import TreeStyle, RectFace, faces, Tree, NodeStyle, TextFace
from ete3 import NCBITaxa
import ast
from Bio import Phylo
import os
import math
# Make sure to install svglib and rlPyCairo - rlPyCairo is not in the imports but is necessary to make svglib and reportlab work
from PIL import Image, ImageFont, ImageDraw
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPM
import time
import numpy as np
import json
class StatisticalFunctions:
        @staticmethod
        def stfdr(p, nBoot=100):
                # set random state
                np.random.RandomState(0)
                
                d = 0.01
                lam = np.arange(0, max(p) + d/2, d)
                m = len(p)
                
                ii = np.argsort(p)
                p = np.sort(p)
                
                pi_0hat = StatisticalFunctions.pi0computer(p,lam,m)
                rem = np.isnan(pi_0hat)
                lam = lam[~rem]
                pi_0hat = pi_0hat[~rem]
                
                mpi_0hat = np.min(pi_0hat)
                pboot = StatisticalFunctions.bootsamp(p, m * nBoot).reshape(m, nBoot)
                
                mse = np.zeros_like(lam)
                for j in range(nBoot):
                        pi_0hatboot = StatisticalFunctions.pi0computer(pboot[:, j], lam, m)
                        mse += (pi_0hatboot - mpi_0hat) ** 2
                
                lm = StatisticalFunctions.LocalMaxMin(mse)
                lmse = lm * (np.max(mse) - mse)
                min_lmse_idx = np.where(lmse == np.min(lmse))[0]
                pi_0 = np.min(pi_0hatboot[min_lmse_idx])
                
                q = StatisticalFunctions.compute_q(p, pi_0, m)
                q[ii] = q
                return q, pi_0
        @staticmethod
        def pi0computer(p, lam, m):
                pi_0hat = np.zeros_like(lam)
                
                for k in range(len(lam)):
                        pi_0hat[k] = np.sum(p > lam[k]) / (m * (1 - lam[k]))
                
                return pi_0hat
        @staticmethod
        def bootsamp(x, ns):
                x = np.array(x).flatten()
                n = len(x)
                s = np.random.rand(ns)
                n_idx = np.digitize(s, np.linspace(0, 1, n+1)) - 1
                s = x[n_idx]
                return s
        @staticmethod
        def LocalMaxMin(vec):
                        scores = np.zeros_like(vec, dtype=int)
                        for i in range(len(vec)):
                                if i == 0:
                                # First element
                                        scores[i] = 0 
                                elif i == len(vec) - 1:
                                # Last element
                                        scores[i] = 0 
                                else:
                                # All other elements
                                        scores[i] = 0
                                        if vec[i] > vec[i-1] and vec[i] > vec[i+1]:
                                                scores[i] = 1
                                        elif vec[i] < vec[i-1] and vec[i] < vec[i+1]:
                                                scores[i] = -1
                        return scores
        @staticmethod
        def compute_q(p, pi_0, m):
                sf = np.arange(1, m+1)
                q = pi_0 * ((p * m) / sf)
                for k in range(m-2, -1, -1):
                        if q[k] > q[k+1]:
                                q[k] = q[k+1]
                return q
class TreeMaker:
        def __init__(self, items_to_find, directorypath, cnodesdir, treetitle):
            #Initialise inputs
            with open(items_to_find) as sf:
                    text = sf.readlines()
                    newls = []
                    for i in text:
                            newls.append(i.rstrip("\n"))
            self.items_to_find = newls
            self.directorypath = directorypath
            self.cnodesdir = cnodesdir
            self.treetitle = treetitle
        def listmaker(self, listtobeprocessed, allitems):
                #Lists all the descendants of items to be found (the txt file which is the user input) and adds them to allitems
                ncbi = NCBITaxa()
                keys = []
                for itf in listtobeprocessed:
                        for i in ncbi.get_lineage(itf):
                                keys.append(str(i))
                allitems = dict.fromkeys(keys, 0)
                newproc = []
                for ltbp in listtobeprocessed:
                        newid = ncbi.get_lineage(ltbp)[::-1][0]
                        newproc.append(str(newid))
                for np in newproc:
                        if np in allitems.keys():
                                allitems[np] +=1
                return allitems

        # JMP's implementation of the ST-FDR (q-value)
        def colourselecter(self, colourdict):
                allitems = self.listmaker(self.items_to_find, {})
                ncbi = NCBITaxa()
                counts = list(allitems.values())
                print(counts)
                counts = np.array(counts).T
                taxa = list(allitems.keys())
                taxrank = []
                for i in taxa: 
                        rank = ncbi.get_rank([i])[int(i)]
                        taxrank.append(rank)
                
                print(taxrank)
                im = np.zeros(len(taxrank), dtype=int)
                mtap = np.zeros_like(im, dtype=float)  # Adjusted for possible shape differences
                mtac = np.zeros_like(im, dtype=float) 
                rnkdict = {}
                txr = list(set(taxrank))
                print(txr)
                for rnk in txr:
                        im = (np.array(taxrank) == rnk).astype(int)
                        counts_im = counts * im  # Element-wise multiplication
                        n = np.sum(counts_im)
                        k = im.sum()
                        if k == 0:
                                continue  # Skip if no elements match the current rank

                        nsamp = 99999  # number of resamplings
                        R = np.zeros((nsamp, k))  # empirical estimate

                        # Generate random samples more efficiently
                        rand_choices = np.random.randint(0, k, (nsamp, int(np.ceil(n))))
                        np.add.at(R, (np.arange(nsamp)[:, None], rand_choices), 1)
                        R /= R.sum(axis=1, keepdims=True)  # Normalize rows
                        R *= n

                        r0 = counts[im == 1]
                        pp = (np.count_nonzero(R > r0, axis=0) + 1) / (nsamp + 1)

                        mtap[im == 1] = pp  # multiple testing adjusted p-values
                        qq = StatisticalFunctions.stfdr(pp)
                        mtac[im == 1] = qq[0]  # multiple testing adjusted q-values
                mtac2 = -np.log10(mtac)
                weighteddict2 = {k:v for k,v in zip(list(taxa), mtac2)}
                print(weighteddict2)
                total = max(mtac2)
                viridis = ['#fde725','#f8e621','#f1e51d','#ece51b','#e5e419','#dfe318','#d8e219','#d0e11c','#cae11f','#c2df23','#bddf26','#b5de2b','#addc30','#a8db34','#a0da39','#9bd93c','#93d741','#8ed645','#86d549','#7fd34e','#7ad151','#73d056','#6ece58','#67cc5c','#60ca60','#5cc863','#56c667','#52c569','#4cc26c','#48c16e','#42be71','#3dbc74','#3aba76','#35b779','#32b67a','#2eb37c','#2ab07f','#28ae80','#25ac82','#24aa83','#22a785','#20a486','#1fa287','#1fa088','#1f9e89','#1e9b8a','#1f998a','#1f968b','#20938c','#20928c','#218f8d','#228d8d','#238a8d','#24878e','#25858e','#26828e','#26818e','#277e8e','#287c8e','#29798e','#2a768e','#2b748e','#2c718e','#2d708e','#2e6d8e','#306a8e','#31688e','#32658e','#33638d','#34608d','#365d8d','#375b8d','#38588c','#39558c','#3b528b','#3c508b','#3d4d8a','#3e4989','#3f4788','#414487','#424186','#433e85','#443a83','#453882','#463480','#46327e','#472e7c','#472c7a','#482878','#482475','#482173','#481d6f','#481b6d','#481769','#471365','#471063','#460b5e','#46085c','#450457','#440154']
                reverseviridis = viridis[::-1]
                #Attributes a number to each colour in the viridis scale for accessing later. 
                colourscaledict = {}
                for i in range(0, len(reverseviridis)):
                        colourscaledict[str(reverseviridis[i])]= i
                tblabelled = []
                # Deciding which taxa should be annotated with names (and added to the tblabelled list)
                for key in weighteddict2.keys():
                        rankdict = ncbi.get_rank(key) 
                        if 'subspecies' in rankdict.values() or 'species' in rankdict.values() or 'genus' in rankdict.values() or 'family' in rankdict.values():

                                lineage = ncbi.get_lineage(key)
                                LEN = len(lineage)
                                nottblabelled = False
                                parent = lineage[LEN-1]
                                parentweight = weighteddict2[str(parent)]
                                descendants = ncbi.get_descendant_taxa(key)
                                descendantweight = 0
                                for i in descendants:
                                        if i in weighteddict2.keys():
                                                descendantweight += weighteddict2[str(i)]
                                        else:
                                                continue
                                if 'subspecies' in rankdict.values() and weighteddict2[key]*100/total>20:
                                        if weighteddict2[key]>parentweight and weighteddict2[key]>descendantweight:
                                                tblabelled.append(key)
                                elif 'species' in rankdict.values() and weighteddict2[key]*100/total >60 :
                                        if weighteddict2[key]>parentweight and weighteddict2[key]>descendantweight:
                                                tblabelled.append(key)
                                elif 'genus' in rankdict.values() and weighteddict2[key]*100/total>80:
                                        if weighteddict2[key]>parentweight or weighteddict2[key]>descendantweight:
                                                tblabelled.append(key)
                                elif 'family' in rankdict.values() and weighteddict2[key]*100/total>100:
                                        if weighteddict2[key]>parentweight and weighteddict2[key]>descendantweight:
                                                tblabelled.append(key)
                        
                          
                        
                for key,value in weighteddict2.items():
                        if value == 0:
                                colourdict[key] = '#440154'
                        else:       
                                placeindex = (value / total) * 100
                                placeindex = math.ceil(placeindex)
                                colourdict[key] = reverseviridis[placeindex - 1]
                penultimate = []
                for tbl in tblabelled:
                        lineagetbl = ncbi.get_lineage(tbl)
                        penultimate.append(lineagetbl[len(lineagetbl)-2])
                parents = list(set(penultimate))
                siblings = {}
                for p in parents:
                        sib = {}
                        for tbl in tblabelled:
                                lineagetbl = ncbi.get_lineage(tbl)
                                parent = (lineagetbl[len(lineagetbl)-2])
                                if parent == p:
                                        sib[tbl] = weighteddict2[tbl]
                        siblings[p] = sib
                tblabelled2 = []
                for key, value in siblings.items():
                        tblabelled2.append(max(value, key = value.get))
                tblabelled2 = list(set(tblabelled2))
                #Save important lists and dictionaries to txtfiles folder
                os.mkdir(self.directorypath + "/" + self.treetitle)
                os.mkdir(self.directorypath + "/" + self.treetitle + "/txtfiles")
                time.sleep(2)
                with open(self.directorypath + "/" + self.treetitle + '/txtfiles/colourdict.txt', 'w', encoding = 'utf-8') as f:
                        f.write(str(colourdict))
                with open(self.directorypath + "/" + self.treetitle + '/txtfiles/total.txt', 'w', encoding = 'utf-8') as g:
                        g.write(str(total))
                with open(self.directorypath + "/" + self.treetitle + '/txtfiles/tblabelled2.txt', 'w', encoding = 'utf-8') as h:
                        h.write(str(tblabelled2))
                with open(self.directorypath + "/" + self.treetitle + '/txtfiles/weighteddict2.txt', 'w', encoding = 'utf-8') as i:
                        i.write(str(weighteddict2))
                #Pause function to prevent the errno 13 error
                time.sleep(2)
                with open(self.directorypath + "/" + self.treetitle + '/txtfiles/colourscaledict.txt', 'w', encoding = 'utf-8') as j:
                        j.write(str(colourscaledict))
                time.sleep(2)
                
                                        

        def layoutfunc(self, node):
                  ncbi = NCBITaxa()
                  rankdict = ncbi.get_rank([node.name])

                  node.complete_branch_lines_when_necessary = False
                  node.optimal_scale_level = "full"
                  node.guiding_lines_type = 0
                  node.extra_branch_line_type = 0
                  with open(self.directorypath + "/" + self.treetitle + '/txtfiles/colourdict.txt') as f:
                                colourdict = f.read()
                  with open(self.directorypath + "/" + self.treetitle + '/txtfiles/tblabelled2.txt') as g:
                                tblabelled = g.read()
                  with open(self.directorypath + "/" + self.treetitle + '/txtfiles/colourscaledict.txt', 'r') as j:
                                indicator = j.read()             
                  colourdict = ast.literal_eval(colourdict)
                  tblabelled = ast.literal_eval(tblabelled)
                  indicator = ast.literal_eval(indicator)
                  tblabellednames = ncbi.get_taxid_translator(tblabelled)
                  
                  
                  node.img_style["hz_line_type"] = 0
                  
                  if node.get_children() == [] or node.name not in colourdict.keys():
                          node.img_style["hz_line_color"] = "#ffffff"
                  nohorline = False
                  

                  if node.name in colourdict.keys():
                        
                        if indicator[colourdict[str(node.name)]] >= 60:
                                amplifier = 15
                        elif indicator[colourdict[str(node.name)]] <=59 and indicator[colourdict[str(node.name)]] >=20:
                                amplifier = 10
                        elif indicator[colourdict[str(node.name)]] <=19 and indicator[colourdict[str(node.name)]] >=3:
                                amplifier = 3
                        else:
                                amplifier = indicator[colourdict[str(node.name)]]
                        node.img_style["fgcolor"] = colourdict[str(node.name)]
                        node.img_style["size"] =0
                        node.img_style["vt_line_color"] = colourdict[str(node.name)]
                        node.img_style["hz_line_color"] = colourdict[str(node.name)]
                        node.img_style["vt_line_width"] = 0
                        node.img_style["hz_line_width"] = 0

                        if node.get_children == []:
                                for i in node.get_children():
                                        if i in colourdict.keys():
                                                continue
                                        else:
                                                nohorline == True
                        
                        if node.name in tblabelled:
                                        if 'species' or 'subspecies' or 'genus' in rankdict.values():
                                                        node.img_style["hz_line_type"] = 0
                                        if amplifier == 15:
                                                textsize = 10
                                                node.img_style["hz_line_type"] = 1
                                                faces.add_face_to_node(TextFace(tblabellednames[int(node.name)], ftype = 'arial', fsize = textsize, fgcolor = "000000", penwidth=0, fstyle= 'normal', tight_text = False, bold = False), node, column = 1, position = "float")
                                        elif amplifier == 10:
                                                textsize = 9
                                                node.img_style["hz_line_type"] = 1
                                                faces.add_face_to_node(TextFace(tblabellednames[int(node.name)], ftype = 'arial', fsize = textsize, fgcolor = "000000", penwidth=0, fstyle= 'normal', tight_text = False, bold = False), node, column = 3, position = "float")
                                        elif amplifier == 3:
                                                textsize = 8
                                                node.img_style["hz_line_type"] = 1
                                                faces.add_face_to_node(TextFace(tblabellednames[int(node.name)], ftype = 'arial', fsize = textsize, fgcolor = "000000", penwidth=0, fstyle= 'italic', tight_text = False, bold = False), node, column = 3, position = "float")
                                        else:
                                                textsize = 0 
                  else: 
                        
                        node.img_style["size"] = 0
                        node.img_style["fgcolor"] = "#ffffff"
                        node.img_style["bgcolor"] = "#ffffff"
                        node.support = 0
                        node.distance = 0
                        node.img_style["vt_line_color"] = "#ffffff"
                        node.img_style["hz_line_color"] = "#ffffff"
                        
                        if node.get_children() != []:
                                nohorline == True
                                        
                  if nohorline == True:
                              node.img_style["hz_line_color"] = "#ffffff"   


                              # Produce the tree
                            
        def Maker(self):     
                                
                                     with open(self.cnodesdir,encoding = 'utf-8') as fp:
                                                      
                                                      text = fp.readlines()
                                                      topologyfeeder = []
                                                      for i in text:
                                                            if " " in i:
                                                                    
                                                                id = i.split(" ")[0]
                                                            else:
                                                                id = i
                                                            topologyfeeder.append(str(id))
                                                
                                                      ncbi = NCBITaxa()
                                                      
                                                      tree = ncbi.get_topology(topologyfeeder, intermediate_nodes=True)
                                                      tree.annotate_ncbi_taxa()
                                                      ts = TreeStyle()
                                                      ts.layout_fn = self.layoutfunc
                                                      ts.show_leaf_name = False
                                                      ts.mode = "c"
                                                      ts.root_opening_factor = 0
                                                      ts.arc_start = 0 # 0 degrees = 3 o'clock
                                                      ts.arc_span = 360
                                                      self.colourselecter({})
                                                      tree.show(tree_style=ts)
                                                      os.mkdir(self.directorypath + "/" + self.treetitle + "/trees")
                                                      tree.write(format = 0, outfile = self.directorypath + "/" + self.treetitle + "/trees/new_tree.nwk")
                                                      tree.render(self.directorypath + "/" + self.treetitle + "/trees/img_tree.svg", w= 3600, units = 'px', tree_style = ts)
                                                      Phylo.convert(self.directorypath + "/" + self.treetitle + "/trees/new_tree.nwk", "newick", self.directorypath + "/" + self.treetitle + "/trees/new_tree.xml", "nexml")
                                                      img = svg2rlg(self.directorypath + "/" + self.treetitle + "/trees/img_tree.svg")
                                                      renderPM.drawToFile(img, self.directorypath + "/" + self.treetitle + "/trees/tree.png", fmt = "PNG")

                                                      filename = self.directorypath + "/" + self.treetitle + "/trees/tree.png"
                                                      with Image.open(filename) as img: 
                                                                width, height = img.size
                                                                img = img.resize((width * 2, height * 2 ))
                                                                img.save(self.directorypath + "/" + self.treetitle + "/trees/tree.png")
                                                                img2 = Image.open('colourbar.png')
                                                                img.paste(img2, (10, 10))
                                                                img.save(self.directorypath + "/" + self.treetitle + "/trees/tree+bar.png")
                                                      with open(self.directorypath + "/" + self.treetitle + '/txtfiles/total.txt') as j:
                                                                        total = j.read()
                                                      img = Image.open(self.directorypath + "/" + self.treetitle + "/trees/tree+bar.png")
                                                      draw = ImageDraw.Draw(img)
                                                      font = ImageFont.truetype("arial", 50)
                                                      font2 = ImageFont.truetype("arial", 70)
                                                      draw.text((600, 150), "Weighted total=  " + str(total), (0, 0, 0), font = font)
                                                      
                                                      img.save(self.directorypath + "/" + self.treetitle  + "/trees/tree+bar.png")
                                                      draw.text((2500, 150), "Phylogenetic Tree based on the named entities in the folder: "  + self.treetitle, (0, 0, 0), font = font2)
                                                      img.save(self.directorypath + "/" + self.treetitle  + "/trees/tree+bar.png")